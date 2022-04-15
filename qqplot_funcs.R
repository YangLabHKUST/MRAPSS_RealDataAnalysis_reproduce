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



qqplot3_sep <- function(p1, p2, p3,
                        label1="Supported",
                        label2="Unknown",
                        label3 ="g3",
                        max=5,strong=NULL){
  
  ci = 0.95
  df <- data.frame(pval = c(sort(p1),
                            sort(p2),
                            sort(p3)),
                   Group = c(rep(label1, length(p1)),
                             rep(label2, length(p2)),
                             rep(label3, length(p3))))
  
  df$pval[-log10(df$pval)>=max] = 10^(-1*max)
  
  label = c(label1, label2, label3)[which.max(c(length(p1), length(p2), length(p3)))]
  
  df <- data.frame(Group = df$Group,
                   observed = -log10(df$pval),
                   expected = c(-log10(ppoints(length(p1))),-log10(ppoints(length(p2))),-log10(ppoints(length(p3)))),
                   clower   = c(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:length(p1), shape2 = (length(p1)):1)),
                                -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:length(p2), shape2 = (length(p2)):1)),
                                -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:length(p3), shape2 = (length(p3)):1))),
                   cupper   = c(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:length(p1), shape2 = (length(p1)):1)),
                                -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:length(p2), shape2 = (length(p2)):1)),
                                -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:length(p3), shape2 = (length(p3)):1))))
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  plot = ggplot(df,aes(expected, observed, color=Group, shape=Group)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(data = subset(df, Group== label), aes(expected, cupper), linetype = 2, color="#999999") +
    geom_line(data = subset(df, Group== label), aes(expected, clower), linetype = 2, color="#999999") +
    scale_shape_manual(name="Trait pair",
                       breaks = c(label1,label2,label3),
                       values = c(16,17,18),
                       labels = c("All", 
                                  expression(paste(r[g]," insignificant")),
                                  expression(paste("Both ", r[g]," and ", c[12]," insignificant"))))+
    scale_color_manual(name="Trait pair",
                       breaks = c(label1,label2,label3),
                       values = c(npg[1],npg[4],npg[3]),
                       labels = c("All", 
                                  expression(paste(r[g]," insignificant")),
                                  expression(paste("Both ",r[g]," and ", c[12]," insignificant"))),
    ) +
    xlab(log10Pe) +
    ylab(log10Po) +
    xlim(0,max)+
    scale_y_continuous(breaks=c(0,1,5), labels=c(0,1,5), limits = c(0,5)) +
    theme_pubr() +
    theme(legend.position = "top",
          legend.text = element_text(hjust = 0, size=20),
          legend.title = element_blank(),
          legend.key.size = unit(2,"line"),
          axis.title = element_text(size = 20, color = "black"),
          axis.text  = element_text(size = 15,color = "black"),
          plot.title = element_text(hjust=0.5, size = 20, face = "bold"))
  
  if(!is.null(strong)) plot= plot + geom_point(data=subset(df, Group==strong),size=1,shape=16)
  
  return(plot)
}

qqplot4 <- function(dat.group,
                    comp1=NULL, comp2=NULL,
                    comp3=NULL, comp4=NULL,
                    comp.labels=NULL,max=45,strong=NULL){
  ci = 0.95
  n  <- nrow(dat.group)
  df <- data.frame(pval = c(sort(dat.group[, comp1]),
                            sort(dat.group[, comp2]),
                            sort(dat.group[, comp3]),
                            sort(dat.group[, comp4])),
                   Method = c(rep(comp1, n),
                              rep(comp2, n),
                              rep(comp3, n),
                              rep(comp4, n)))
  
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
                       levels = c(comp1, comp2,comp3,comp4)
    )
  }else{
    df$Method = factor(df$Method,
                       levels = c(comp1, comp2,comp3,comp4),
                       labels = comp.labels
    )
  }
  
  
  plot = ggplot(df,aes(expected, observed, color=Method)) +
    geom_point(size = 3, shape = 16) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, color="#999999") +
    geom_line(aes(expected, clower), linetype = 2, color="#999999") +
    xlab(log10Pe) +
    ylab(log10Po) +
    ylim(0,max) +
    guides(color = guide_legend(override.aes = list(size = 2)))+
    scale_color_manual(values=c("#F8766D","#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_pubr() +
    theme(legend.position = c(1,0.025),
          legend.justification = c(1.1,0),
          legend.text = element_text(hjust = 0, size=15),
          legend.title = element_text(hjust = 0, size=15),
          legend.key.size = unit(1,"line"),
          axis.title = element_text(size = 20, color = "black"),
          axis.text  = element_text(size = 15,color = "black"),
          plot.title = element_text(hjust = 0.5, size = 20, vjust=2)
    )
  
  if(!is.null(strong)) plot= plot + geom_point(data=subset(df, Method==strong),size=1,shape=16)
  return(plot)
}

