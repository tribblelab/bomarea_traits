library(RColorBrewer)

thinLogs = function(log.file.path,out.file.path=NA,keep.or.discard,n.species) {
  log = read.table(log.file.path,header=T)
  keep = c(T,T,T,T,rep(keep.or.discard,n.species),rep(T,(length(keep.or.discard)^2)))
  if(is.na(out.file.path)) {
    out.file.path = paste(strsplit(log.file.path,".log")[[1]][1],"_thinned.log",sep="")
  }
  write.table(log[,keep],out.file.path,row.names=F,col.names=T,quote=F,sep="\t")
}

extractPrecision = function(log.file.path,out.file.path=NA) {
  log = read.table(log.file.path,header=T)
  keep = c(1,grep("precision",names(log)))
  if(is.na(out.file.path)) {
    out.file.path = paste(strsplit(log.file.path,".log")[[1]][1],"_precision.log",sep="")
  }
  write.table(log[,keep],out.file.path,row.names=F,col.names=T,quote=F,sep="\t")
}

extractVarCovar = function(log.file.path,out.file.path=NA,trait.names=NA) {
  # recover()
  header <- readLines(con=log.file.path,n=5)
  header <- header[!grepl("#",header)]
  header <- header[1]
  header <- strsplit(header,"\t")[[1]]
  keep <- grep("precision",header)
  ntraits <- sqrt(length(keep))
  if(is.na(trait.names[1])) {
    trait.names <- c(1:ntraits)
  } else if (length(trait.names) != ntraits) {
  	stop(paste("Please provide ",ntraits," trait names\n"))
  }
  count = 0
   col.names <- character(ntraits^2)
   for (i in 1:ntraits) {
      for (j in 1:ntraits) {
        count = count + 1
        col.names[count] = paste(trait.names[i],"_with_",trait.names[j],sep="")
      }
    }
  log = read.table(log.file.path,header=T)
  # keep = grep("precision",names(log))
  # ntraits = sqrt(length(keep))
  precision = log[,keep]
  nobs = length(precision[,1])
  var.covar = matrix(ncol=length(keep),nrow=nobs)
  
  for (count in 1:nobs) {
    var_covar = solve(matrix(as.numeric(precision[count,]),ncol=ntraits,nrow=ntraits,byrow=TRUE))
    vec = unlist(var_covar)
    var.covar[count,] = vec
  }
  
  if(is.na(out.file.path)) {
    out.file.path = paste(strsplit(log.file.path,".log")[[1]][1],"_var_covar.log",sep="")
  }
    
  var.covar = cbind(log[,1],var.covar)
  
  write.table(var.covar,out.file.path,col.names=c("state", col.names),row.names=F,quote=F,sep="\t")
}


extractCorrelations = function(log.file.path,out.file.path=NA,trait.names=NA) {
  log = read.table(log.file.path,header=T)
  keep = grep("precision",names(log))
  ntraits = sqrt(length(keep))
  precision = log[,keep]
  nobs = length(precision[,1])
  correlations = matrix(ncol=length(keep),nrow=nobs)
  
  for (count in 1:nobs) {
  	var_covar = solve(matrix(as.numeric(precision[count,]),ncol=ntraits,nrow=ntraits,byrow=TRUE))
    vec = NULL
    for (i in 1:ntraits) {
   	  for (j in 1:ntraits) {
  	    vec = c(vec, var_covar[i,j]/sqrt(var_covar[i,i]*var_covar[j,j]) )
  	  }
    }
    correlations[count,] = vec
  }
    
  if(is.na(out.file.path)) {
    out.file.path = paste(strsplit(log.file.path,".log")[[1]][1],"_correlations.log",sep="")
  }
  
  if(is.na(trait.names[1])) {
  	count = 0
  	correlation.names = NULL
  	for (i in 1:ntraits) {
  		for (j in 1:ntraits) {
  			count = count + 1
  			correlation.names[count] = paste(i,"_with_",j,sep="")
  		}
  	}
  } else {
  	if ( length(trait.names) != length(keep) && length(trait.names) != sqrt(length(keep)) ) {
  		return(paste("There are ",ntraits," traits, so there are ",ntraits**2,"correlations,. Please provide either ",ntraits," trait names or ",ntraits**2," correlation names in input \"trait.names\""))
  	}
  	if ( length(trait.names) == sqrt(length(keep)) ) {
	  	count = 0
	  	correlation.names = NULL
	  	for (i in 1:ntraits) {
	  		for (j in 1:ntraits) {
	  			count = count + 1
	  			correlation.names[count] = paste(trait.names[i],"_with_",trait.names[j],sep="")
	  		}
	  	}
	  }
  }
  
  correlations = cbind(log[,1],correlations)
  
  write.table(correlations,out.file.path,col.names=c("state",correlation.names),row.names=F,quote=F,sep="\t")
}


plotCorrelationDistributions = function(correlations,labels,lwd=1,cex.labels=1,omi=c(1,1,1,1)) {
  
  ntraits = sqrt(ncol(correlations))
  
  # now we will get column indices of each unique trait comparison
  col.index = rep(1:ntraits,ncol(correlations)/ntraits)
  row.index = NULL
  for (i in 1:ntraits) {
    row.index = c(row.index,rep(i,ntraits))
  }
  
  keepers = 1:ncol(correlations)
  keepers = keepers[which(row.index != col.index)]
  keepers = keepers[which(row.index[keepers] < col.index[keepers])]
 
  corr.means = apply(correlations[,keepers],2,mean)
  
  # assign colors 
  neg.corr.colors = brewer.pal(9,"Blues")
  pos.corr.colors = brewer.pal(9,"Reds")
  
  pos = which(corr.means > 0)
  neg = which(corr.means < 0)
  
  colors = NULL
  fill.colors = NULL
 
  colors[pos] = pos.corr.colors[cut(corr.means[pos],seq(0,1,1/6), labels = F)+3]
  fill.colors[pos] = pos.corr.colors[cut(corr.means[pos],seq(0,1,1/6), labels = F)]
  colors[neg] = neg.corr.colors[cut(-corr.means[neg],seq(0,1,1/6), labels = F)+3]
  fill.colors[neg] = neg.corr.colors[cut(-corr.means[neg],seq(0,1,1/6), labels = F)]
  
  for (i in 1:choose(ntraits,2)) { # now we grey-out the non-significant correlations
    samples <- correlations[,keepers[i]]
    q <- quantile(samples,prob=c(0.025,0.975))
    
    if(q[1] > 0 && q[2] > 0){ # positive correlation
    } else if (q[1] < 0 && q[2] < 0){ # negative correlation
    } else { # no correlation
      colors[i] <- 'black'
      fill.colors[i] <- '#7f7f7f50'
    }
  }
  
  layout.mat <- matrix(0,nrow=ntraits,ncol=ntraits)
  layout.mat[lower.tri(layout.mat)] <- 1:(choose(ntraits,2))
  layout.mat = t(layout.mat)
  layout.mat[!upper.tri(layout.mat)] <- (choose(ntraits,2)+1):ncol(correlations)

  layout(layout.mat)
  par(mar=c(0,0,0,0),omi=omi,xpd=T)
 
   for (count in 1:choose(ntraits,2)) {
      samples <- correlations[,keepers[count]]
      dens <- density(samples,from=-1,to=1)
      plot(dens,xlim=c(-1,1),xaxt='n',yaxt='n',zero.line=FALSE,main=NA,xlab=NA,ylab=NA,lwd=lwd,col=colors[count])
      polygon(dens,border=NA,col=fill.colors[count])
      abline(v=0,lty=2)
    }
  
  count = 0
  for (j in 1:ntraits) {
    for (i in j:ntraits) {
      plot(NULL,xlim=c(-1,1),ylim=c(-1,1),xaxt='n',yaxt='n',main=NA,xlab=NA,ylab=NA,bty="n")
      if (i != j) {
        count = count + 1
        
        polygon(x=c(-1,-1,1,1),y=c(-1,1,1,-1),col=fill.colors[count],border=NA)
        text(labels=round(corr.means[count],3),x=0,y=0)
      }
      if (i == ntraits) {
         mtext(labels[j],side=1,line=1,cex=cex.labels)
      }
      if (j == 1) {
        mtext(labels[i],side=2,las=1,line=2,cex=cex.labels)
      }
    }
  }
}

# library(RColorBrewer)
# external <- read.table("~/git_repos/Andy/spermetazoa/data_files/traces/state-specific_correlations/ginsi_external_combined_correlations.log",row.names=1,header=T,stringsAsFactors=F)
# internal <- read.table("~/git_repos/Andy/spermetazoa/data_files/traces/state-specific_correlations/ginsi_internal_combined_correlations.log",row.names=1,header=T,stringsAsFactors=F)
# labels = c("internal fertilization","feeding larvae","egg size","sperm-body length","sperm-body width")
# 
# blue = brewer.pal(9,"Set1")[2]
# orange = brewer.pal(9,"Set1")[5]
# 
# 
# compareCorrelationDistributions(external[],internal[],labels=rev(labels[-1]),background.upper=paste(blue,"30",sep=""),background.lower=paste(orange,"30",sep=""))

compareCorrelationDistributions = function(correlations.upper,correlations.lower,labels,lwd=1,cex.labels=1,omi=c(1,1,1,1),background.upper=NULL,background.lower=NULL) {
  recover()
  ntraits = sqrt(ncol(correlations.upper))
  
  # now we will get column indices of each unique trait comparison
  col.index = rep(1:ntraits,ncol(correlations.upper)/ntraits)
  row.index = NULL
  for (i in 1:ntraits) {
    row.index = c(row.index,rep(i,ntraits))
  }
  
  keepers = 1:ncol(correlations.upper)
  keepers = keepers[which(row.index != col.index)]
  keepers = keepers[which(row.index[keepers] < col.index[keepers])]
  
  corr.means.upper = apply(correlations.upper[,keepers],2,mean)
  corr.means.lower = apply(correlations.lower[,keepers],2,mean)
  
  # assign colors for upper diagonal
  neg.corr.colors = brewer.pal(9,"Blues")
  pos.corr.colors = brewer.pal(9,"Reds")
  
  pos.upper = which(corr.means.upper > 0)
  neg.upper = which(corr.means.upper < 0)
  
  colors.upper = NULL
  fill.colors.upper = NULL
  
  colors.upper[pos.upper] = pos.corr.colors[cut(corr.means.upper[pos.upper],seq(0,1,1/6), labels = F)+3]
  fill.colors.upper[pos.upper] = pos.corr.colors[cut(corr.means.upper[pos.upper],seq(0,1,1/6), labels = F)]
  colors.upper[neg.upper] = neg.corr.colors[cut(-corr.means.upper[neg.upper],seq(0,1,1/6), labels = F)+3]
  fill.colors.upper[neg.upper] = neg.corr.colors[cut(-corr.means.upper[neg.upper],seq(0,1,1/6), labels = F)]
  
  for (i in 1:choose(ntraits,2)) { # now we grey-out the non-significant correlations
    samples <- correlations.upper[,keepers[i]]
    q <- quantile(samples,prob=c(0.025,0.975))
    
    if(q[1] > 0 && q[2] > 0){ # positive correlation
    } else if (q[1] < 0 && q[2] < 0){ # negative correlation
    } else { # no correlation
      colors.upper[i] <- 'black'
      fill.colors.upper[i] <- '#7f7f7f50'
    }
  }
  
  # assign colors for lower diagonal
  neg.corr.colors = brewer.pal(9,"Blues")
  pos.corr.colors = brewer.pal(9,"Reds")
  
  pos.lower = which(corr.means.lower > 0)
  neg.lower = which(corr.means.lower < 0)
  
  colors.lower = NULL
  fill.colors.lower = NULL
  
  colors.lower[pos.lower] = pos.corr.colors[cut(corr.means.lower[pos.lower],seq(0,1,1/6), labels = F)+3]
  fill.colors.lower[pos.lower] = pos.corr.colors[cut(corr.means.lower[pos.lower],seq(0,1,1/6), labels = F)]
  colors.lower[neg.lower] = neg.corr.colors[cut(-corr.means.lower[neg.lower],seq(0,1,1/6), labels = F)+3]
  fill.colors.lower[neg.lower] = neg.corr.colors[cut(-corr.means.lower[neg.lower],seq(0,1,1/6), labels = F)]
  
  for (i in 1:choose(ntraits,2)) { # now we grey-out the non-significant correlations
    samples <- correlations.lower[,keepers[i]]
    q <- quantile(samples,prob=c(0.025,0.975))
    
    if(q[1] > 0 && q[2] > 0){ # positive correlation
    } else if (q[1] < 0 && q[2] < 0){ # negative correlation
    } else { # no correlation
      colors.lower[i] <- 'black'
      fill.colors.lower[i] <- '#bfbfbf'
    }
  }
  
  layout.mat <- matrix(0,nrow=ntraits,ncol=ntraits)
  layout.mat[lower.tri(layout.mat)] <- 1:(choose(ntraits,2))
  layout.mat <- t(layout.mat)
  layout.mat[lower.tri(layout.mat)] <- (choose(ntraits,2)) + 1:(choose(ntraits,2))
  diag(layout.mat) <- (2*(choose(ntraits,2))+1):ntraits**2
  
  layout(layout.mat)
  par(mar=c(0,0,0,0),omi=omi,xpd=T)
  
  for (count in 1:choose(ntraits,2)) {
    samples <- correlations.upper[,keepers[count]]
    dens <- density(samples,from=-1,to=1)
    plot(dens,xlim=c(-1,1),xaxt='n',yaxt='n',main=NA,xlab=NA,ylab=NA,type="n")
    y.usr <- c(par("usr")[3],par("usr")[4],par("usr")[4],par("usr")[3])
    x.usr <- c(par("usr")[1],par("usr")[1],par("usr")[2],par("usr")[2])
    polygon(x=x.usr,y=y.usr,border=NA,lty=0,col=background.upper)
    lines(dens,xlim=c(-1,1),xaxt='n',yaxt='n',main=NA,xlab=NA,ylab=NA,lwd=lwd,col=colors.upper[count])
    polygon(dens,border=NA,col=fill.colors.upper[count])
    abline(v=0,lty=2)
  }
  
  for (count in 1:choose(ntraits,2)) {
    samples <- correlations.lower[,keepers[count]]
    dens <- density(samples,from=-1,to=1)
    plot(dens,xlim=c(-1,1),xaxt='n',yaxt='n',main=NA,xlab=NA,ylab=NA,type="n")
    y.usr <- c(par("usr")[3],par("usr")[4],par("usr")[4],par("usr")[3])
    x.usr <- c(par("usr")[1],par("usr")[1],par("usr")[2],par("usr")[2])
    polygon(x=x.usr,y=y.usr,border=NA,lty=0,col=background.lower)
    lines(dens,xlim=c(-1,1),xaxt='n',yaxt='n',main=NA,xlab=NA,ylab=NA,lwd=lwd,col=colors.lower[count])
    # plot(dens,xlim=c(-1,1),xaxt='n',yaxt='n',zero.line=FALSE,main=NA,xlab=NA,ylab=NA,lwd=lwd,col=colors.lower[count])
    polygon(dens,border=NA,col=fill.colors.lower[count])
    abline(v=0,lty=2)
  }
}



plotLatentLiabilityVarianceComparison <- function(variances,discrete.or.continuous,lwd=1,cex.labels=1,omi=c(0,0,0,0),mai=c(0.25,0.25,0.25,0.25),trait.names=NULL,trait.colors=NULL) {
  # recover()
  if (class(variances) != "list") {
    if (is.null(trait.names[1])) {
      trait.names <- names(variances.all)
    }
    variances.all <- as.list(variances)
  } else {
    if (length(trait.names) != length(variances)) {
      stop('If input is of format "list", argument "trait.names" must be provided, as a vector containing names for each trait')
    }
    variances.all <- variances
  }
  if (length(discrete.or.continuous) != length(variances)) {
    stop('You must use argument "discrete.or.continuous" to specify whether each column/list element in "variances" refers to a discrete or continuous trait')
  }
  discrete.or.continuous <- tolower(discrete.or.continuous)
  n.continuous <- length(grep("c",discrete.or.continuous))
  n.discrete <- length(grep("d",discrete.or.continuous))
  if (n.continuous > 0) {
    if (n.discrete > 0) {
      n.panels <- 2
    } else {
      n.panels <- 1
    }
  } else {
    n.panels <- 1
  }
  if (n.panels == 2) {
    layout(rbind(c(1,2)))
  }
  par(omi=omi,mai=mai,xpd=T)
  if (n.continuous > 0) {
    continuous.legend.names <- trait.names[which(discrete.or.continuous == "c")]
    variances <- variances.all[which(discrete.or.continuous == "c")]
    if (is.null(trait.colors[1])) {
      if (n.continuous > 2) {
        auto.colors <- brewer.pal(n.continuous,"Set2")
        colors <- auto.colors
      } else if (n.continuous == 2) {
        auto.colors <- brewer.pal(9,"Set2")[c(2,4)]
        colors <- auto.colors
      } else if (n.continuous == 1) {
        stop("There is only 1 continuous trait. There is no value in plotting a single trait for comparison")
      }
    } else {
      auto.colors <- NULL
      if (n.continuous == 1) {
        stop("There is only 1 continuous trait. There is no value in plotting a single trait for comparison")
      } else {
        colors <- trait.colors[which(discrete.or.continuous == "c")]
      }
    }
    min.x <- floor(min(unlist(variances)))
    max.x <- ceiling(max(unlist(variances)))
    max.dens <- numeric(n.continuous)
    for (i in 1:n.continuous) {
      samples <- variances[[i]]
      dens <- density(samples)
      max.dens[i] <- max(dens$y)
    }
    plot(NULL,NULL,xlim=c(min.x,max.x),ylim=c(0,max(max.dens)),main="",xlab="",ylab="")
    for (i in 1:n.continuous) {
      samples <- variances[[i]]
      dens <- density(samples)
      lines(dens,xlim=c(min.x,max.x),ylim=c(0,max(max.dens)),xaxt='n',yaxt='n',main=NA,xlab=NA,ylab=NA,lwd=lwd,col=colors[i])
      polygon(dens,border=NA,col=paste(colors[i],"95",sep=""))
    }
    legend(x="topright",y=NULL,legend=continuous.legend.names,fill=colors,bty="n",border=NA)
    mtext("density",side=2,line=1.25,outer=T)
  }
  if (n.discrete > 0) {
    discrete.legend.names <- trait.names[which(discrete.or.continuous == "d")]
    variances <- variances.all[which(discrete.or.continuous == "d")]
    if (is.null(trait.colors[1])) {
      if (n.discrete > 2) {
        auto.colors <- brewer.pal(n.discrete,"Set2")
        colors <- auto.colors
      } else if (n.discrete == 2) {
        auto.colors <- brewer.pal(9,"Set2")[c(2,4)]
        colors <- auto.colors
      } else if (n.discrete == 1) {
        stop("There is only 1 discrete trait. There is no value in plotting a single trait for comparison")
      }
    } else {
      auto.colors <- NULL
      if (n.discrete == 1) {
        stop("There is only 1 discrete trait. There is no value in plotting a single trait for comparison")
      } else {
        colors <- trait.colors[which(discrete.or.continuous == "d")]
      }
    }
    min.x <- floor(min(unlist(variances)))
    max.x <- ceiling(max(unlist(variances)))
    max.dens <- numeric(n.discrete)
    for (i in 1:n.discrete) {
      samples <- variances[[i]]
      dens <- density(samples)
      max.dens[i] <- max(dens$y)
    }
    plot(NULL,NULL,xlim=c(min.x,max.x),ylim=c(0,max(max.dens)),main="",xlab="",ylab="")
    for (i in 1:n.discrete) {
      samples <- variances[[i]]
      dens <- density(samples)
      lines(dens,xlim=c(min.x,max.x),ylim=c(0,max(max.dens)),xaxt='n',yaxt='n',main=NA,xlab=NA,ylab=NA,lwd=lwd,col=colors[i])
      polygon(dens,border=NA,col=paste(colors[i],"95",sep=""))
    }
    legend(x="topright",y=NULL,legend=discrete.legend.names,fill=colors,bty="n",border=NA)
    if (n.continuous == 0) {
      mtext("density",side=2,line=1.25,outer=T)
    }
  }
}
