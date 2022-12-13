MP_gLasso <- function(cv_object, group, lambda.type = "min", sort.type = "mean", max.shown = 20, intercept=TRUE) {
  
  ## extracting value ##
  ## extracting value - coefficients and lambda ##
  
  ifelse(lambda.type == "min", coef <- coef(cv_object$gglasso.fit, s=cv_object$lambda.min),
         coef <- coef(cv_object$gglasso.fit, s=cv_object$lambda.1se))
  ifelse(lambda.type == "min", lambda <- round(cv_object$lambda.min, 4),
         lambda <- round(cv_object$lambda.1se, 4))
  
  if (intercept) {
    coef_name <- dimnames(coef)[[1]][-1]
    coef <- coef[-1]
  } else {coef_name <- dimnames(coef)[[1]]}
  
  if (mode(group) == "numeric"){
    group <- as.character(group)
  }
  
  idx <- (coef != 0)
  idx2 <- group %in% group[idx]
  
  beta <- coef[idx2]
  beta_name <- coef_name[idx2]
  
  group_name <- group[idx2]
  group_size <- as.integer(table(group_name))
  
  gname <- unique(group_name)
  gid <- as.numeric(factor(group_name))
  
  len <- length(gname)
  
  if(len == 0){
    stop('No variable with non-zero coefficient')
  }
  if(len> max.shown){ 
    if(sort.type == 'mean'){
        avgbeta_abs <- tapply(abs(beta), group_name, mean)  
        oo <- order(avgbeta_abs, decreasing = TRUE)
        gname = rownames(avgbeta_abs[oo][1:max.shown])
        idx2 = group %in% gname
        beta <- coef[idx2]
        beta_name <- coef_name[idx2]
      
        group_name <- group[idx2]
        group_size <- as.integer(table(group_name))
      
        gname <- unique(group_name)
        gid <- as.numeric(factor(group_name))
      
        len <- length(gname)
    }else{
      max.beta_abs <- tapply(abs(beta), group_name, max)  
      oo <- order(max.beta_abs, decreasing = TRUE)
      gname = rownames(max.beta_abs[oo][1:max.shown])
      idx2 = group %in% gname
      beta <- coef[idx2]
      beta_name <- coef_name[idx2]
      
      group_name <- group[idx2]
      group_size <- as.integer(table(group_name))
      
      gname <- unique(group_name)
      gid <- as.numeric(factor(group_name))
      
      len <- length(gname)
    }
  }

  if (sort.type == "mean") {
    avgbeta_abs <- tapply(abs(beta), group_name, mean)  
    oo <- order(avgbeta_abs, decreasing = TRUE)
    scale <- tapply(abs(beta), group_name, mean) / tapply(abs(beta), group_name, max)
    
    data <- data.frame(beta=beta, group_name=factor(gid))
    data$sign <- ifelse(data$beta > 0, "beta \nwith positive sign", "beta \nwith negative sign")
    
    p2 <- data %>%
      group_by(group_name) %>%
      summarise(avg_beta = mean(abs(beta)), gsize=n()) %>%
      mutate(group_name = fct_reorder(group_name, avg_beta, .desc=TRUE),
             group_size = fct_reorder(factor(gsize), gsize, .desc=FALSE)) %>%
      ggplot() +
      geom_bar(aes(x=group_name, y=avg_beta, fill=group_size), col="black",
               stat="identity", width = 1) +
      theme_light() +
      scale_fill_grey(start=1, end=0.5) +
      scale_x_discrete(breaks=factor(unique(gid)), labels=gname) +
      coord_polar()
  } else {
    max.beta_abs <- tapply(abs(beta), group_name, max)  
    oo <- order(max.beta_abs, decreasing = TRUE)
    scale <- rep(1, length(unique(group_name)))
    
    data <- data.frame(beta=beta, group_name=factor(gid))
    data$sign <- ifelse(data$beta > 0, "beta \nwith positive sign", "beta \nwith negative sign")
    
    p2 <- data %>%
      group_by(group_name) %>%
      summarise(max_beta = max(abs(beta)), gsize=n()) %>%
      mutate(group_name = fct_reorder(group_name, max_beta, .desc=TRUE),
             group_size = fct_reorder(factor(gsize), gsize, .desc=FALSE)) %>%
      ggplot() +
      geom_bar(aes(x=group_name, y=max_beta, fill=group_size), col="black",
               stat="identity", width = 1) +
      theme_light() +
      scale_fill_grey(start=1, end=0.5) +
      scale_x_discrete(breaks=factor(unique(gid)), labels=gname) +
      coord_polar()
  }
  
  position <- list(length = len)
  for (i in 1:len) {
    position[[i]] <- jitter(rep(i, group_size[oo][i]))
  }
  
  beta_scale <- list(length = len)
  for (i in 1:len) {
    beta_scale[[i]] <- split(abs(beta), group_name)[[i]] * scale[i]
  }
  
  beta_sign <- list(length = len)
  for (i in 1:len) {
    beta_sign[[i]] <- split(data$sign, group_name)[[i]]
  }
  
  beta_origin <- list(length = len)
  for (i in 1:len) {
    beta_origin[[i]] <- split(beta, group_name)[[i]]
  }
  
  indexx <- list(length = len)
  for (i in 1:len){
    indexx[[i]] <- split(seq(sum(group_size)), group_name)[[i]]
  }
  
  b_name <- list(length = len)
  for (i in 1:len) {
    b_name[[i]] <- split(beta_name, group_name)[[i]]
  }
  
  sample_points <- list(length= len)
  for (i in 1:len) {
    sample_points[[i]] <- data.frame(position=position[[i]], beta_scale=beta_scale[[oo[i]]],
                                     beta_sign=beta_sign[[oo[i]]], index=indexx[[oo[i]]], coef=beta_origin[[oo[i]]],
                                     name=b_name[[oo[i]]])
  }
  
  for (i in 1:len) {
    sample_points[[i]]$tt <- paste0("Name : ", sample_points[[i]]$name,
                                    "\nCoefficient : ", round(sample_points[[i]]$coef,4))
  }
  
  p3 <- p2
  for(i in 1:len){
    loop_input <- paste0('geom_point_interactive(data=sample_points[[i]], aes(x=position, y=beta_scale,
                       shape=beta_sign, colour=beta_sign, tooltip=tt, data_id=index))')
    p3 <- p3 + eval(parse(text = loop_input))
  }
  
  p4 <- p3 +
    scale_shape_manual(name = 'sign of beta', values=c(1, 17)) +
    scale_colour_manual(name = 'sign of beta', values = c('black', 'black', 'black')) +
    labs(title = "Group Lasso",
         subtitle = paste0("Lambda : ", lambda, ", Lambda type : ",lambda.type,"\nSort type : ", sort.type),
         x="", y="", fill="Group size") +
    theme(plot.title = element_text(hjust = 0.5,size=13, face="bold"), 
          legend.key = element_rect(colour = "black"),
          plot.subtitle = element_text(size=11, face="bold")) +
    guides(fill = guide_legend(ncol=2)) +
    scale_y_continuous(labels = NULL)
  
  girafe(code = print(p4), width_svg = 7, height_svg = 6)
  
  
}

