comple <- function(allele){

  ifelse(allele == "A","T", ifelse(allele == "T","A", ifelse(allele == "G","C", ifelse(allele == "C","G", allele)) ))

}

harmonise_ref <- function(dat1, dat2){

  SNPorder = match(dat2$SNP, dat1$SNP)
  dat1 = dat1[SNPorder,]
  dat1$A1 <-  as.character(dat1$A1)
  dat1$A2 <-  as.character(dat1$A2)

  dat2$A1 <- as.character(dat2$A1)
  dat2$A2 <-  as.character(dat2$A2)

  flip.index = which((dat1$A1 == dat2$A2 & dat2$A1 == dat1$A2) |
                       (dat1$A1 ==comple(dat2$A2) & dat2$A1 == comple(dat1$A2)))

  dat1$A1 <-  dat2$A1
  dat1$A2  <- dat2$A2

  dat1[flip.index ,"b.exp"] = - dat1[flip.index ,"b.exp"]
  dat1[flip.index ,"b.out"] = - dat1[flip.index ,"b.out"]

  return(dat1)
}
