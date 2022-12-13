MP_Lasso <- function(cv_object, lambda.type = "min", max.shown = 20, intercept=TRUE) {
  
  ## extracting value ##
  ## extracting value - coefficients and lambda ##
  
  ifelse(lambda.type == "min", coef <- coef(cv_object$glmnet.fit, s=cv_object$lambda.min),
         coef <- coef(cv_object$glmnet.fit, s=cv_object$lambda.1se))
  ifelse(lambda.type == "min", lambda <- round(cv_object$lambda.min, 4),
         lambda <- round(cv_object$lambda.1se, 4))
  
  if (intercept) {
    coef_name <- dimnames(coef)[[1]][-1]
    coef <- coef[-1]
  } else {coef_name <- dimnames(coef)[[1]]}
  
  idx <- (coef != 0)

  beta <- coef[idx]
  beta_name <- coef_name[idx]
  
  len = length(beta)
  if(len == 0){
    stop('No variable with non-zero coefficient')
  }
  if(len > max.shown){
    idx = order(abs(beta), decreasing = TRUE)[1:max.shown]
    beta = beta[idx]
    beta_name = beta_name[idx]
  }
  data = data.frame(name = beta_name, beta = abs(beta))
  data$sign = ifelse(beta > 0 , 'positive', 'negative')
  data$name = factor(beta_name, levels = unique(beta_name))
  data$tt = paste0('Name: ', data$name, '\nCoefficient : ', round(data$beta, 4))
  p2 = ggplot(data) + 
    geom_bar_interactive(aes(x=name, y=beta, fill = sign, tooltip = tt), color = 'black', stat="identity", width = 1) +
    theme_light() +
    scale_fill_grey(start=1, end=0.5) +
    scale_x_discrete(breaks=factor(unique(beta_name)), labels=unique(beta_name))  + 
    coord_polar()
  
  p3 <- p2 +
    labs(title = "Lasso",
         subtitle = paste0("Lambda : ", lambda, ", Lambda type : ",lambda.type),
         x="", y="", fill="sign of beta") +
    theme(plot.title = element_text(hjust = 0.5,size=13, face="bold"), 
          legend.key = element_rect(colour = "black"),
          plot.subtitle = element_text(size=11, face="bold")) +
    guides(fill = guide_legend(ncol=2)) +
    scale_y_continuous(labels = NULL)
  
  girafe(code = print(p3), width_svg = 7, height_svg = 6)
}

