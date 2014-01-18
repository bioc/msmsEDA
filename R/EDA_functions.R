###################################
###     PACKAGE    msmsEDA      ###
###      SpC EDA PIPE-LINE      ###
###        by J. Gregori        ###
###    v. 12th January, 2014    ###
###################################

###  Pre-process msms counts matrix by replacing NAs by 0s, 
###    and further removing all-zero and '-R' pprotein rows
pp.msms.data <- function(msnset)
{
  ###  Set 0 for NA
  counts <- exprs(msnset)
  if (sum(is.na(counts))) 
  { counts[is.na(counts)] <- 0
    exprs(msnset) <- counts
  }
  ###  Remove all zero rows and '-R' proteins
  fl1 <- apply(counts, 1, function(x) sum(x) > 0)
  pnms <- rownames(counts)
  fl2 <- substring(pnms, nchar(pnms) - 1) == "-R"
  flags <- (!fl2 & fl1)
  msnset <- msnset[flags, ]
  MSnbase:::logging(msnset,"Applied pp.msms.data preprocessing")
}


###  Gene names and accession table  
gene.table <- function(Accession,Protein,patt="GN=[A-Z0-9_]*",off=3)
{ ich <- regexpr(patt,Protein)
  nch <- attr(ich,"match.length")
  nafl <- ich<0
  nch[nafl] <- 2; ich[nafl] <- 0
  gn.tbl <- substr(Protein,ich+off,ich+nch-1)
  gn.tbl[nafl] <- "NA"
  names(gn.tbl) <- Accession
  gn.tbl
}


###  SpC statistics by sample  
count.stats <- function(msnset)
{ msms.counts <- exprs(msnset)
  #  the five numbers for SpC by sample 
  fvnm <- t(apply(msms.counts,2,fivenum))
  colnames(fvnm) <- c("min","lwh","med","hgh","max") 
  #  identified proteins by sample
  pps <- apply(msms.counts,2,function(x) sum(x>0))
  #  total counts by sample
  cps <- apply(msms.counts,2,sum)
  tfvnm <- cbind(proteins=pps,counts=cps,fvnm)
  tfvnm
}


###  Samples PCA map
counts.pca <- function(msnset,facs=NULL,do.plot=TRUE,snms=NULL,wait=TRUE)
{ #  PCA decomposition
  msms.counts <- exprs(msnset)
  smpl.pca <- prcomp(t(msms.counts))
  #  PC variances
  pc.vars <- summary(smpl.pca)$importance
  if(!do.plot)
    invisible( list(pca=smpl.pca,pc.vars=pc.vars) )
  #  Plots
  cls <- c(brewer.pal(8,"Dark2"),"red","navy","black")
  if(is.null(snms))
    snms <- paste("X",1:ncol(msms.counts),sep="")
  xlb <- paste("PC1  ",round(summary(smpl.pca)$importance[2,1]*100,2),
               "% var")
  ylb <- paste("PC2  ",round(summary(smpl.pca)$importance[2,2]*100,2),
               "% var")
  if(is.null(facs)) facs <- pData(msnset)
  if(!is.data.frame(facs)) facs <- data.frame(facs)
  for(i in 1:ncol(facs))
  { if(i>1 & interactive() & wait)
      readline(prompt = "Pause. Press <Enter> to continue...")
    gtit <- paste("PCA - ",colnames(facs)[i],sep="")
    eqscplot(smpl.pca$x[,1],smpl.pca$x[,2],type="n",xlab= xlb,ylab=ylb)
    text(smpl.pca$x[,1],smpl.pca$x[,2],snms,cex=0.7,
         font=2,col=cls[as.integer(facs[,i])])
    abline(h=0,v=0,lty=4,col="gray")
    title(main=gtit,line=1)
  }
  invisible( list(pca=smpl.pca,pc.vars=pc.vars) )
}
	  

# Highlight the edges bringing to a set of labels
edgeCol <- function(dend, keys, fgr="red", lwd=1, ...) 
{ myattr <- attributes(dend)
  if(is.leaf(dend)) 
  { # Si es una fulla amb etiqueta en keys
    if(length(which(keys==myattr$label))==1)
      attr(dend,"edgePar") <- c(myattr$edgePar, 
           list(col=fgr,lwd=lwd))
  } else {
    lbl <- labels(dend)
    # Si es una branca amb totes les fulles en keys
    if( length(intersect(keys,lbl))==length(lbl) )
      attr(dend,"edgePar") <- c(myattr$edgePar,
           list(col=fgr,lwd=lwd))
  }
  return(dend)
}


###  Samples dendrograms
counts.hc <- function(msnset,do.plot=TRUE,facs=NULL,wait=TRUE)
{ msms.counts <- exprs(msnset)  
  hc <- hclust(dist(t(msms.counts)))
  if(!do.plot)
    invisible(hc)
	
  smpl.nms <- colnames(msms.counts)
  par(cex.main=1,cex.lab=0.8,cex.axis=0.8)
  omar <- par(mar=c(4.5,2,2.5,5))
  cls <- c(brewer.pal(8,"Dark2"),"red","navy","black")
  if(is.null(facs)) facs <- pData(msnset)
  if(!is.data.frame(facs)) facs <- data.frame(facs)
  for(i in 1:ncol(facs))
  { if(i>1 & interactive() & wait)
      readline(prompt = "Pause. Press <Enter> to continue...")
    gtit <- paste("HC - ",colnames(facs)[i],sep="")
    dend_colo1 <- as.dendrogram(hc)
    for(j in 1:nlevels(facs[,i]))
       dend_colo1 <- dendrapply(dend_colo1,edgeCol, fgr=cls[j], 
               lwd=2,keys=smpl.nms[facs[,i]==levels(factor(facs[,i]))[j]]) 
    plot(dend_colo1,nodePar=list(pch=NA,lab.cex=0.7),horiz=TRUE)
    title(main=gtit,line=1,cex=1)
  } 
  par(mar=omar)
  invisible(hc)
}


###  Normalize sample counts by given divisors
norm.counts <- function(msnset,div)
{ msms.counts <- exprs(msnset)
  exprs(msnset) <- sweep(msms.counts,MARGIN=2,STATS=div,FUN="/")
  MSnbase:::logging(msnset,"Applied norm.counts normalization ")
}


###  Samples heatmap 
counts.heatmap <- function(msnset,etit=NULL,fac=NULL,to.pdf=FALSE)
{ msms.norm <- exprs(msnset)
  if(to.pdf)
    pdf(file=paste(etit,"-HeatMap.pdf",sep=""),width=6,height=6)
    
  if( is.null(etit) ) etit <- "SpC"
  if( !is.null(fac) )
  { cls <- c("red","blue")[as.integer(fac)]
    hm <- heatmap(t(scale(t(msms.norm))),col=greenred(255),labRow=NA,
                  cexCol=0.7,ColSideColors=cls)
   } else
      hm <- heatmap(t(scale(t(msms.norm))),col=greenred(255),labRow=NA,
                    cexCol=0.7)
  #  whole size heatmap (3mm by protein in height)
  if(to.pdf)
  { dev.off()
    h <- nrow(msms.norm)/(2.54/0.3)
    pdf(file=paste(etit,"-FullHeatMap.pdf",sep=""),
        width=7,height=h)
    if( !is.null(fac) )
    { cls <- c("red","blue")[as.integer(fac)]
      hm <- heatmap.2(t(scale(t(msms.norm))),col=greenred(255),trace="none",
                    key=FALSE,cexRow=0.6,cexCol=0.7,margins=c(5,6),
				    dendrogram="both",lhei=c(2,h-3),ColSideColors=cls)
    } else
      heatmap.2(t(scale(t(msms.norm))),col=greenred(255),trace="none",
                key=FALSE,cexRow=0.6,cexCol=0.7,margins=c(5,6),
                dendrogram="both",lhei=c(2,h-3))
    dev.off()
  }
}


###  Residual variance of a one factor fit.
residual.var <- function(msms,gpf)
{ ss <- function(x)
          var(x)*(length(x)-1)
  apply(msms,1,function(x)
         sum(tapply(x,gpf,ss))/(ncol(msms)-nlevels(gpf)) )
}

###  One factor mean of level means
class.means <- function(msms,gpf)
{ 
  apply(msms,1,function(x) mean(tapply(x,gpf,mean)))
}

###  Compute and plot distribution of dispersion coefficients
###    by considering one factor at a time. 
disp.estimates <-
       function(msnset,facs=NULL,do.plot=TRUE,etit=NULL,to.pdf=FALSE,wait=TRUE)
{ msms.norm <- exprs(msnset)
  if(is.null(facs)) facs <- pData(msnset)
  if(!is.data.frame(facs)) facs <- data.frame(facs)
  if(is.null(etit)) etit <- "MSMS"
  if(do.plot & to.pdf)
  { pdf(file=paste(etit,"-DispPlots.pdf",sep=""),width=6,height=10,
        paper="a4")
    par(cex.main=1,cex.lab=0.8,cex.axis=0.8)
    par(mfrow=c(2,1))
  }
  p <- c(0.25,0.5,0.75,0.9,0.95,0.99,1)
  res <- matrix(0,nrow=ncol(facs),ncol=length(p))
  colnames(res) <- p
  rownames(res) <- colnames(facs)
  for(i in 1:ncol(facs))
  { pmn <- class.means(msms.norm,facs[,i])
    pvar <- residual.var(msms.norm,facs[,i])
    pdisp <- pvar/pmn
    pdisp[is.na(pdisp)] <- 0
    res[i,] <- quantile(pdisp,prob=p)
    if(do.plot)
    { if(!to.pdf & i>1 & interactive() & wait)
        readline(prompt = "Pause. Press <Enter> to continue...")
      gtit <- paste("Residual dispersion factor ",colnames(facs)[i],sep="")
      plot(density(pdisp),xlim=c(0,quantile(pdisp,0.975)),main="")
      abline(v=1,lty=4,col="gray")
      title(main=gtit,line=1,cex=1)
      if(!to.pdf & interactive() & wait)
        readline(prompt = "Pause. Press <Enter> to continue...")
      eqscplot(log10(pmn),log10(pvar),pch="+",cex=0.6,xlab="Log SpC mean",
               ylab="Log SpC variance")
      abline(a=0,b=1,lty=4,col="blue")
      title(main=gtit,line=1,cex=1)

    }
  }
  if(to.pdf) dev.off()
  invisible(res)
}


###  Flag features of low signal and lowest dispersion
filter.flags <- function(data,minSpC=2,frac.out=0.4)
{ data <- data[rowMeans(data)>=minSpC,]
  cdisp <- apply(data,1,var)/rowMeans(data)
  thr <- quantile(cdisp,probs=frac.out)
  cdisp>=thr
}


###  SpC barplots by sample
spc.barplots <- function(msms.counts,fact=NULL,...)
{ n <- ncol(msms.counts)
  if(is.null(fact))
    fact <- rep(1,n)
  if(is.data.frame(fact) | is.matrix (fact)) fact <- fact[,1]
  fact <- as.factor(fact)   
  pal <- brewer.pal(8,"Dark2")
  cls <- pal[as.integer(fact)]
  tcnts <- apply(msms.counts,2,sum)
  medt <- median(tcnts)
  div <- tcnts/medt
  names(div) <- colnames(msms.counts)
  omar <- par(mar=c(7,3,2,2))
  if("col" %in% names(list(...)))
  { barplot(div,las=2,...)
  } else {
    barplot(div,las=2,col=cls,...)
  }
  abline(h=1,lty=2,col="blue")
  par(mar=omar)
}


###  SpC boxplots by sample
spc.boxplots <- function(msms.counts,fact=NULL,minSpC=2,...)
{ n <- ncol(msms.counts)
  if(is.null(fact))
     fact <- rep(1,n)
  if(is.data.frame(fact) | is.matrix (fact)) fact <- fact[,1]
  fact <- as.factor(fact)   
  
  ###  Boxplots
  mx <- 0
  slst <- vector("list",n)
  names(slst) <- colnames(msms.counts)
  fn <- function(x) log2(x+0.1)
  for(i in 1:n)
  { v <- msms.counts[,i]
    slst[[i]] <- fn(v[v>=minSpC])
    mx <- max(mx,slst[[i]])
  }
  omar <- par(mar=c(6.5,4,4,2)+0.1)
  if("col" %in% names(list(...)))
  { pal <- list(...)$col
    cls <- pal[as.integer(fact)]
    boxplot(slst,ylim=c(1,mx*1.15),xaxt="n",ylab="log2 SpC",...)
  } else {
    pal <- brewer.pal(8,"Dark2")
    cls <- pal[as.integer(fact)]
    boxplot(slst,ylim=c(1,mx*1.15),xaxt="n",ylab="log2 SpC",col=cls,...)
  }
  axis(side=1,at=1:n,las=2,cex.axis=0.8,labels=names(slst))
  if(minSpC)
  { tt <- paste("min SpC >=",minSpC)
    legend("topright",fill=pal,legend=levels(fact),cex=0.8,
           title=tt)
  } else {           
    legend("topright",fill=pal,legend=levels(fact),cex=0.8)
  }
  par(mar=omar)
}


###  SpC density plots by sample
spc.densityplots <- function(msms.counts,fact=NULL,minSpC=2,...)
{ 
  n <- ncol(msms.counts)
  if(is.null(fact))
     fact <- rep(1,n)
  if(is.data.frame(fact) | is.matrix (fact)) fact <- fact[,1]
  fact <- as.factor(fact)   
  pal <- brewer.pal(8,"Dark2")
  cls <- pal[as.integer(as.factor(fact))]
  ###  Density plots
  mx <- 0
  n <- ncol(msms.counts)
  slst <- vector("list",n)
  fn <- function(x) log2(x+0.1)
  for(i in 1:n)
  { v <- msms.counts[,i]
    slst[[i]] <- fn(v[v>=minSpC])
    mx <- max(mx,density(slst[[i]])$y)
  }
  omar <- par(mar=c(5,4,3,2)+0.1)
  if("main" %in% names(list(...)))
  { plot(density(slst[[i]]),ylim=c(0,mx),col=cls[1],lwd=2,
           xlab="log2 SpC",...)
  } else {
    plot(density(slst[[i]]),ylim=c(0,mx),col=cls[1],lwd=2,
         xlab="log2 SpC",main="",...)
  }
  if(minSpC>0)
    abline(v=log2(minSpC),lty=4,col="gray")
  for(i in 1:n)
    lines(density(slst[[i]]),col=cls[i],lwd=2)
  if(minSpC)
  { tt <- paste("min SpC >=",minSpC)
    legend("topright",fill=pal,legend=levels(fact),cex=0.8,
           title=tt)
  } else {           
    legend("topright",fill=pal,legend=levels(fact),cex=0.8)
  }
  par(mar=omar)
}


###  Scatterplot means vs means given a two level factor
###  Two SpC transformations are possible, either log2 or sqrt.
spc.scatterplot <- function(msms.counts,treat,trans="log2",minSpC=2,minLFC=1,...)
{
  # signal and effect filter
  treat <- as.factor(treat)   
  mspc <- t(apply(msms.counts,1,function(x) tapply(x,treat,mean)))
  mx <- apply(mspc[,1:2],1,max)
  mn <- apply(mspc[,1:2],1,min)
  fl <- apply(mspc[,1:2],1,function(x) max(x)>=minSpC) & abs(log2(mx/mn))>=minLFC
  # scatterplot on the untransformed counts
  if(trans=="none")
  { omar <- par(mar=c(5,5,3,2)+0.1)
    plot(mspc[,1:2],pch="+",cex=0.7,asp=1,
       xlab=paste("mean SpC(",colnames(sqrtm)[1],")",sep=""),
       ylab=paste("mean SpC(",colnames(sqrtm)[2],")",sep=""),...)
    points(mspc[fl,1:2],pch="+",cex=0.7,col="red")
    abline(a=0,b=1,lty=4,col="gray")
    abline(a=0,b=2,lty=4,col="navy")
    abline(a=0,b=0.5,lty=4,col="navy")
    par(mar=omar)
  }
  # scatterplot on the sqrt transformed counts
  if(trans=="sqrt")
  { sqrtm <- sqrt(mspc)
    omar <- par(mar=c(5,5,3,2)+0.1)
    plot(sqrtm[,1:2],pch="+",cex=0.7,asp=1,
       xlab=substitute(sqrt(m*e*a*n~ ~SpC(x)), list(x = colnames(sqrtm)[1])),
       ylab=substitute(sqrt(m*e*a*n~ ~SpC(x)), list(x = colnames(sqrtm)[2])),...)
    points(sqrtm[fl,1:2],pch="+",cex=0.7,col="red")
    abline(a=0,b=1,lty=4,col="gray")
    abline(a=0,b=sqrt(2),lty=4,col="navy")
    abline(a=0,b=1/sqrt(2),lty=4,col="navy")
    par(mar=omar)
  }
  # scatterplot on the log2 transformed counts
  off <- 0.01
  if(trans=="log2")
  { logm <- log2(mspc+off)
    mxx <- max(logm[,1])
    mxy <- max(logm[,2])
    mx <- max(mxx,mxy)
    omar <- par(mar=c(5,5,3,2)+0.1)
    plot(logm[,1:2],pch="+",cex=0.7,asp=1,xlim=c(-3,mx),ylim=c(-3,mx),
       xlab=substitute(log[2](m*e*a*n~ ~SpC(x)), list(x = colnames(logm)[1])),
       ylab=substitute(log[2](m*e*a*n~ ~SpC(x)), list(x = colnames(logm)[2])),...)
    points(logm[fl,1:2],pch="+",cex=0.7,col="red")
    abline(a=0,b=1,lty=4,col="gray")
    abline(a=log2(2),b=1,lty=4,col="navy")
    abline(a=-log2(2),b=1,lty=4,col="navy")
    par(mar=omar)
  }
}


###  Batch effects correction
###    Fit a model with a batch factor
###    Remove batch effect. 
###    If half==TRUE, half correction to each batch
###    If sqrt.trans is TRUE transform previously the SpC by sqrt
###     and back transform after the correction.
batch.neutralize <- function(dat,fbatch,half=TRUE,sqrt.trans=TRUE)
{ Z <- dat
  #  Sqrt transformation
  if(sqrt.trans) Z <- sqrt(dat)
  #  The model matrix  with a single factor (batch)
  if(!half)
    ocont <- options(contrasts=c("contr.treatment","contr.poly"))
  if(half)
    ocont <- options(contrasts=c("contr.sum","contr.poly"))
  n <- ncol(dat)
  X <- model.matrix(~fbatch)
  options(contrasts=ocont$contrasts)
  #  Fit the linear model
  XinvXX <- X %*% solve(t(X) %*% X)
  B <- Z %*% XinvXX
  #  The correction 
  dat <- Z - B[,-1,drop=FALSE] %*% t(X[,-1,drop=FALSE])
  #  Back to the original scale
  if(sqrt.trans) dat <- dat^2
  dat
}
