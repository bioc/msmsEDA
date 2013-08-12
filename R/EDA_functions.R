###################################
###     PACKAGE    msmsEDA      ###
###  MS/MS EDA PIPE-LINE v 1.0  ###
###       by J. Gregori         ###
###      24th May, 2013         ###
###################################

###  Add a text line to the MSnSet log and check for validity
add.log <- function(msnset,tproc)
{ msnset@processingData@processing <-
      c(msnset@processingData@processing,paste(tproc,": ",date(),sep=""))
  if (validObject(msnset))
    return(msnset)
}


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
  add.log(msnset,"Applied pp.msms.data preprocessing")
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
counts.pca <- function(msnset,facs=NULL,do.plot=TRUE,snms=NULL)
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
  { if(i>1 & interactive())
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
counts.hc <- function(msnset,do.plot=TRUE,facs=NULL)
{	msms.counts <- exprs(msnset)  
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
  { if(i>1 & interactive())
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
  add.log(msnset,"Applied norm.counts normalization ")
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
function(msnset,facs=NULL,do.plot=TRUE,etit=NULL,to.pdf=FALSE)
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
    res[i,] <- quantile(pdisp,prob=p)
    if(do.plot)
    { if(!to.pdf & i>1 & interactive())
        readline(prompt = "Pause. Press <Enter> to continue...")
      gtit <- paste("Residual dispersion factor ",colnames(facs)[i],sep="")
      plot(density(pdisp),xlim=c(0,quantile(pdisp,0.975)),main="")
      abline(v=1,lty=4,col="gray")
      title(main=gtit,line=1,cex=1)
      if(!to.pdf & interactive())
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
