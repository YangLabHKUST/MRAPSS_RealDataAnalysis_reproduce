qqplot1 <- function(dat.group, comp1, comp.labels=NULL, max=5){
  
  ci = 0.95
  n  <- nrow(dat.group)
  df <- data.frame(pval = c(sort(dat.group[, comp1])),
                   Method = c(rep(comp1, n)))
  
  df$pval[-log10(df$pval)>=max] = 10^(-1*max)
  
  df <- data.frame(Method = df$Method,
                   observed = -log10(df$pval),
                   expected = rep(-log10(ppoints(n)), length(unique(df$Method))),
                   clower   = rep(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))),
                   cupper   = rep(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))))
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  
  
  if(is.null(comp.labels)){
    df$Method = factor(df$Method, levels = comp1)
  }else{
    df$Method = factor(df$Method,
                       levels = comp1,
                       labels = comp.labels)
  }
  
  plot = ggplot(df) +
    geom_point(aes(expected, observed, color=Method), size = 3, shape = 16) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po) +
    ylim(0, max)+
    xlim(0, max)+
    theme_pubr() +
    theme(legend.position = c(1,0.025),
          legend.justification = c(0.95,0),
          legend.text = element_text(hjust = 0, size=20),
          legend.title = element_blank(),
          axis.title = element_text(size = 20, color = "black"),
          axis.text  = element_text(size = 15,color = "black"),
          plot.title = element_text(hjust=-0.1, size = 20, face = "bold"))
  return(plot)
}


qqplot2 <- function(dat.group, comp1, comp2, comp.labels=NULL,max=40, strong=NULL, mark=T){
  ci = 0.95
  n  <- nrow(dat.group)
  df <- data.frame(pval = c(sort(dat.group[, comp1]),
                            sort(dat.group[, comp2])),
                   Method = c(rep(comp1, n),
                              rep(comp2, n)))
  df$pval[-log10(df$pval)>=max] = 10^(-1*max)
  df <- data.frame(Method = df$Method,
                   observed = -log10(df$pval),
                   expected = rep(-log10(ppoints(n)), length(unique(df$Method))),
                   clower   = rep(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))),
                   cupper   = rep(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))))
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  if(is.null(comp.labels)){
    df$Method = factor(df$Method,
                       levels = c(comp1, comp2)
    )
  }else{
    df$Method = factor(df$Method,
                       levels = c(comp1, comp2),
                       labels = comp.labels
    )
  }
  
  plot = ggplot(df, aes(x=expected, y=observed)) +
    geom_point(data = df,
               aes(x=expected, y=observed,color=Method),size = 3, shape = 16) +
    # geom_point(data = subset(df, observed >= -log10(0.05/26/5)),
    #            aes(x=expected, y=observed,color=Method), size = 3, shape = 16) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, color="#999999") +
    geom_line(aes(expected, clower), linetype = 2, color="#999999") +
    xlab(log10Pe) +
    ylab(log10Po) +
    ylim(0, max)+
    xlim(0, max)+
    theme_pubr() +
    theme(legend.position = c(1,0.025),
          legend.justification = c(0.95,0),
          legend.text = element_text(hjust = 0, size=20),
          legend.title = element_blank(),
          axis.title = element_text(size = 20, color = "black"),
          axis.text  = element_text(size = 15,color = "black"),
          plot.title = element_text(hjust=-0.1, size = 20, face = "bold"))
  if(mark) plot = plot +
                  geom_point(data = subset(df, observed > -log10(0.05/26/5)),
                              aes(x=expected, y=observed), color = "darkred", size = 5, shape = 1, stroke=1) +
  
  if(!is.null(strong)) plot= plot + geom_point(data=subset(df, Method==strong),aes(color=Method), size=2,shape=16)
  
  return(plot)
}

