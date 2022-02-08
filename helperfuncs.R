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
  
  out <- list(meta=metadat)
  if ("roi" %in% type) out$roi <- roidat
  if ("stim" %in% type) out$stim <- stimdat
  if ("beh" %in% type) out$beh <- behdat
  
  return(out)
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
      out[[l]] <- rbind(out[[l]], ddat[[l]])
    }
  }
  
  return(out)
}


binnify <- function(roi,stim,wsize,wnum,past=T) {
  slicedat <- roi[stim[,slicetime(soundFrames,wsize,wnum,past),by=.(day,run,stimNum,soundID)],on=.(day,run,t),allow.cartesian=T]
  out <- slicedat[,lapply(.SD,mean),by=.(day,run,stimNum,roi,soundID,winid)]
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


permnull <- function(dat,nsamp=1e3,in_roi=F,tstype="prop_sig") {
  ts <- vector(length=nsamp)
  
  for (i in 1:nsamp) {
    cat(i,"\r")
    
    if (in_roi) {
      dat[,px:=delta*sample(c(1,-1),.N,T)]
    } else {
      dat[,px:=delta*sample(c(1,-1),1),by=stimNum]
    }
    
    ts[i] <- switch(tstype, 
                    pop_mean = dat[,t.test(px)$statistic,by=roi][,t.test(V1)$statistic],
                    prop_sig = dat[,t.test(px)$p.value,by=roi][,mean(V1>0.95)],
                    lrt = logLik(lmer(px ~ 1 + (1|roi) + (1|stimNum),data=dat,REML = F)) - 
                      logLik(lmer(px ~ 0 + (1|stimNum),data=dat,REML = F)))
  }
  
  return(ts)
}

lossybin <- function(ldat,wsize,s=1) {
  bigt <- ldat[,max(t)]
  trim <- ldat[,(bigt-s+1) %% wsize]
  bindat <- ldat[t %in% s:(bigt-trim), .(t=mean(t),lum=mean(lum)), 
       by=.(day,run,roi,cut_interval(t,length=wsize) |> as.numeric())][,-"as.numeric",with=F]
  return(bindat)
}

Wcalc <- function(ldat,wsize=5,s=1,method="ZCA-cor") {
  bdat <- lossybin(ldat,wsize,s)
  ddat <- bdat[,.(t=t[-1],delta=diff(log(lum))),by=.(day,run,roi)] |> dcast(t ~ roi, value.var="delta")
  X <- ddat[,-1,with=F] |> as.matrix()
  mu <- colMeans(X)
  
  if (nrow(X) < ncol(X)) {
    
  }
  
  W <- whiteningMatrix(cov(X),method)
  return(list(mu=mu,W=W))
}


whiten_timepoints <- function(ldat,stim,wsize,s,method="ZCA-cor") {
  W <- Wcalc(ldat,wsize,s,method)
  # if (any(is.na(W$W))) W <- Wcalc(ldat,1,s,method)
  lsdat <- binnify(ldat,stim,wsize,1)
  dsdat <- lsdat[,.(delta=diff(log(lum))),by=.(day,run,stimNum,soundID,roi)]
  return(dsdat[,.(roi,delta,z=as.vector(W$W %*% (delta-W$mu))), by=.(day,run,stimNum,soundID)])
}
