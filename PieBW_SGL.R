##################
## Load library ##
##################
library(gglasso)

library(ggplot2)
library(dplyr)
library(forcats)
library(gridExtra)
library(SGL)
library(pclogit)
library(ggiraph)
library(patchwork)

### example data ##
set.seed(1010)
n <- 100
p <- 200
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
beta <- c(
  rep(3, 4), c(8, -1, 1, 0, 0), rep(-1, 6), c(1, -3, 3, 0, 0),
  rep(0, (p - 20))
)
groups <- c(rep(1,4), rep(2,5), rep(3,6), rep(4,5) ,rep(5:(p / 5), each = 5))
groups_name <- paste0("Group : ", groups)

eps <- rnorm(n, mean = 0, sd = 1)
y <- X %*% beta + eps
pr <- 1 / (1 + exp(-X %*% beta))
y0 <- rbinom(n, 1, pr)

## real data##
############################################

#### Grourp Lasso - regression / binary ####
cv_sgl <- cvSGL(data=list(x=X,y=y), index=groups, type="linear", nfold=3, alpha=0.5, standardize=FALSE)
cv_sgl2 <- cvSGL(data=list(x=X,y=y0), index=groups, type="logit", nfold=3, alpha=0.5, standardize=FALSE)


PieBW_SGL <- function(cv_object=NULL, group = NULL, coef_label = NULL, lambda.type = "min", sort.type = "mean") {
  
  ## extracting value ##
  ## extracting value - coefficients and lambda ##
  
  if (is.null(coef_label)) {
    coef_name <- paste0("V", 1:length(group))
  } else {
    coef_name <- coef_label
  }
  
  df <- apply(cv_object$fit$beta, 2, function(x) sum(x != 0))
  
  ##negative log likelihood
  idx_min <- which(cv_object$lldiff == min(cv_object$lldiff))
  idx_temp <- which((cv_object$lldiff <= min(cv_object$lldiff)+sd(cv_object$lldiff)), 
                    (cv_object$lldiff >= min(cv_object$lldiff)-sd(cv_object$lldiff)))
  idx_1se <- which(df == min(df[idx_temp]))[1]
  
  
  ##coefficients
  ifelse(lambda.type == "min", coef <- cv_object$fit$beta[,idx_min], coef <- cv_object$fit$beta[,idx_1se])
  
  
  ##lambda value##
  ifelse(lambda.type == "min", lambda <- round(cv_object$lambdas[idx_min], 4), 
         lambda <- round(cv_object$lambdas[idx_1se], 4))
  
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
  if (sort.type == "mean") {
    avgbeta_abs <- tapply(abs(beta), group_name, mean)  
    oo <- order(avgbeta_abs, decreasing = TRUE)
    scale <- tapply(abs(beta), group_name, mean) / tapply(abs(beta), group_name, max)
    
    # temp <- rep(1:len, each = as.integer(table(group)))
    # gid <- temp[idx2]
    
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
    labs(title = "Sparse Group Lasso",
         subtitle = paste0("Lambda : ", lambda, ", Lambda type : ",lambda.type,"\nSort type : ", sort.type),
         x="", y="", fill="Group size") +
    theme(plot.title = element_text(hjust = 0.5,size=13, face="bold"), 
          legend.key = element_rect(colour = "black"),
          plot.subtitle = element_text(size=11, face="bold")) +
    guides(fill = guide_legend(ncol=2)) +
    scale_y_continuous(labels = NULL)
  
  girafe(code = print(p4), width_svg = 7, height_svg = 6)
  
  
}


PieBW_SGL(cv_object = cv_sgl, group = groups_name, lambda.type = "min", sort.type = "max")
PieBW_SGL(cv_object = cv_sgl, group = groups_name, lambda.type = "min", sort.type = "mean")
PieBW_SGL(cv_object = cv_sgl2, group = groups_name, lambda.type = "1se", sort.type = "max")
PieBW_SGL(cv_object = cv_sgl2, group = groups_name, lambda.type = "1se", sort.type = "mean")


