dayreader <- function(path, cond="All", type=c("roi","stim","pil","beh")) {
  brbstr <- str_split_fixed(dir(path),"_",4)
  brbrun <- unique(brbstr[,3])
  filebase <- paste(dirname(path) |> basename(), basename(path),sep="_")
  roidat <- data.table()
  metadat <- data.table()
  stimdat <- data.table()
  behdat <- data.table()
  
  for (r in brbrun) {
    brbpre <- file.path(path, paste(filebase,r,sep="_"))
    rmetadat <- fread(paste0(brbpre,"_info.csv"))
    if (!("All" %in% cond) & !(rmetadat$condition %in% cond)) next
    metadat <- rbind(metadat, rmetadat)
    
    if ("roi" %in% type) {
      rroi <- fread(paste0(brbpre,"_Fsoma.csv")) |> t() |> as.data.table()
      rroi[,t:=seq_len(.N)]
      rroi <- melt(rroi, id.vars = "t", variable.name = "roi", value.name = "lum")
      # rroi[,lumz:=scale(lum),by=roi]
      
      if ("pil" %in% type) {
        rpil <- fread(paste0(brbpre,"_Fneuropil.csv")) |> t() |> as.data.table()
        rpil[,t:=1:.N]
        rpil <- melt(rpil, id.vars = "t", variable.name = "roi", value.name = "pil")
        rroi <- rroi[rpil, on=.(t,roi)]
      }
      
      rroi$run <- r
      roidat <- rbind(roidat, rroi)
    }
    
    if ("stim" %in% type) {
      stimtim <- fread(paste0(brbpre,"_stim.csv"))
      stimtim <- cbind(stimNum=seq_len(nrow(stimtim)),stimtim)
      stimtim$run <- r
      
      stimdat <- rbind(stimdat,stimtim)
    }
    
    if ("beh" %in% type) {
      beh <- fread(paste0(brbpre,"_behavior.csv"))
      beh[,t:=seq_len(.N)]
      beh$run <- r
      
      behdat <- rbind(behdat,beh)
    }
    
  }
  
  if (length(metadat)>0) {
    metadat[,run:=str_extract(run,"\\d*$")]
    out <- list(meta=metadat)
    if ("roi" %in% type) out$roi <- roidat
    if ("stim" %in% type) out$stim <- stimdat
    if ("beh" %in% type) out$beh <- behdat
    return(out)
  }
  
}

slicetime <- function(t,wsize,wnum,past=T) {
  trng <- t:(wsize*wnum+t-1)
  winid <- cut_number(trng,n = wnum) |> as.numeric()
  if (past) { 
    winid <- c(sort(-winid),winid-1)
    trng <- c((t-wsize*wnum):(t-1),trng)
  } else {
    winid <- winid - 1
  }
  
  return(data.table(t=trng,winid=winid))
}

brbreader <- function(path, cond="All", type=c("roi","stim","pil","beh")) {
  brbday <- dir(file.path(path))
  out <- list(meta=data.table())
  for (l in type) {
    if (l == "pil") next
    out[[l]] <- data.table()
  }
  
  for (d in brbday) {
    ddat <- dayreader(file.path(path,d), cond, type)
    if (length(ddat[[1]])>0) for (l in names(ddat)) {
      if (l == "pil") next
      ddat[[l]]$day <- d
      out[[l]] <- rbind(out[[l]], ddat[[l]], fill=T)
    }
  }
  
  return(out)
}


binnify <- function(roi,stim,wsize,wnum,past=T) {
  slicedat <- roi[stim[,slicetime(soundFrames,wsize,wnum,past),by=.(day,run,stimNum,soundFrequency,soundID)],on=.(day,run,t),allow.cartesian=T]
  out <- slicedat[,lapply(.SD,mean),by=.(day,run,stimNum,roi,winid)]
  # if (past) out[,postStim:=winid>=0]
  if (out[,any(is.na(roi))]) out <- out[-which(is.na(roi))]
  return(out)
}

whiten <- function(ldat, center=T, scale=T) {
  long <- dcast(ldat,t ~ roi, value.var="lum")
  pcout <- prcomp(long[,-"t",with=F],center=center,scale=scale,retx = T)
  out <- cbind(t=long$t,as.data.table(pcout$x %*% diag(1/pcout$sdev))) |> 
    melt(id.vars = "t", variable.name = "roi", value.name = "lum") |> 
    cbind(run=ldat$run,day=ldat$day)
  
  return(out)
}


# permnull <- function(dat,nsamp=1e3,in_roi=F,tstype="prop_sig") {
# ts <- vector(length=nsamp)
permnull <- function(value,group=NULL) {
  
  if (is.null(group)) {
    pval <- value*sample(c(1,-1),length(value),T)
  } else {
    dat <- data.table(i=1:length(value),value,group)
    dat <- dat[,.(i,value*sample(c(1,-1),1)),by=group]
    setkey(dat,"i")
    pval <- dat$V2
  }
  
  return(pval)
}

lossybin <- function(ldat,wsize,s=1) {
  bigt <- ldat[,max(t)]
  trim <- ldat[,(bigt-s+1) %% wsize]
  bindat <- ldat[t %in% s:(bigt-trim), .(t=mean(t),lum=mean(lum)), 
                 by=.(day,run,roi,cut_interval(t,length=wsize) |> as.numeric())][,-"as.numeric",with=F]
  return(bindat)
}

Wcalc <- function(ldat,wsize=5,s=1,method="ZCA-cor",cov.est="shrinkage") {
  bdat <- lossybin(ldat,wsize,s)
  ddat <- bdat[,.(t=t[-1],delta=diff(log(lum))),by=.(day,run,roi)] |> dcast(t ~ roi, value.var="delta")
  X <- ddat[,-1,with=F] |> as.matrix()
  mu <- colMeans(X)
  
  if (cov.est=="shrinkage") {
    Shat <- cov.shrink(X,lambda.var = 0)
  } else {
    Shat <- cov(X)
  }
  
  W <- whiteningMatrix(Shat,method)
  return(list(mu=mu,W=W))
}

whiten_timepoints <- function(ldat,stim,wsize,s,method="ZCA-cor",w0size=wsize,cov.est="shrinkage") {
  W <- Wcalc(ldat,w0size,s,method,cov.est)
  # if (any(is.na(W$W))) W <- Wcalc(ldat,1,s,method)
  lsdat <- binnify(ldat,stim,wsize,1)
  dsdat <- lsdat[,.(delta=diff(log(lum))),by=.(day,run,stimNum,soundID,roi)]
  return(dsdat[,.(roi,delta,z=as.vector(W$W %*% (delta-W$mu))), by=.(day,run,stimNum,soundID)])
}

pcday <- function(ldat,scale=T,rank=5,renorm=F) {
  X <- dcast(ldat, run + t ~ roi, value.var="lum")
  pca <- prcomp(X[,-(1:2),with=F],center=T,scale=scale,rank=rank)
  Y <- pca$x
  if (renorm) Y <- Y %*% diag(1/pca$sdev[1:rank])
  colnames(Y) <- paste0("PC",1:rank)
  Y <- cbind(X[,.(run,t)],as.data.table(Y))
  return(melt(Y,id.vars=c("run","t"), variable.name="component", value.name="lum"))
}

pcrun <- function(ldat,scale=T,rank=5,renorm=F) {
  X <- dcast(ldat, t ~ roi, value.var="lum")
  pca <- prcomp(X[,-1,with=F],center=T,scale=scale,rank=rank)
  Y <- pca$x
  if (renorm) Y <- Y %*% diag(1/pca$sdev[1:rank])
  colnames(Y) <- paste0("PC",1:rank)
  Y <- cbind(t=X[,t],as.data.table(Y))
  return(melt(Y,id.vars="t", variable.name="component", value.name="lum"))
}