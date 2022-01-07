dayreader <- function(path, cond="All", type=c("roi","stim","pil")) {
  brbstr <- str_split_fixed(dir(path),"_",4)
  brbrun <- unique(brbstr[,3])
  filebase <- paste(dirname(path) |> basename(), basename(path),sep="_")
  roidat <- data.table()
  metadat <- data.table()
  stimdat <- data.table()
  pildat <- data.table()
  
  for (r in brbrun) {
    brbpre <- file.path(path, paste(filebase,r,sep="_"))
    rmetadat <- fread(paste0(brbpre,"_info.csv"))
    if (!("All" %in% cond) & !(rmetadat$condition %in% cond)) next
    metadat <- rbind(metadat, rmetadat)
    
    if ("roi" %in% type) {
      rroi <- fread(paste0(brbpre,"_Fsoma.csv")) |> t() |> as.data.table()
      rroi[,t:=1:.N]
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
    
    
    
  }
  
  out <- list(meta=metadat)
  if ("roi" %in% type) out$roi <- roidat
  if ("stim" %in% type) out$stim <- stimdat
  
  return(out)
}

slicetime <- function(t,wsize,wnum) {
  trng <- t:(wsize*wnum+t-1)
  winid <- cut_number(trng,n = wnum) |> as.numeric()
  winid <- c(sort(-winid),winid-1)
  trng <- c((t-wsize*wnum):(t-1),trng)
  
  return(data.table(t=trng,winid=winid))
}

brbreader <- function(path, cond="All", type=c("roi","stim","pil")) {
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


binnify <- function(roi,stim,wsize,wnum) {
  slicedat <- roi[stim[,slicetime(soundFrames,wsize,wnum),by=.(day,run,stimNum,soundID)],on=.(day,run,t),allow.cartesian=T]
  out <- slicedat[,lapply(.SD,mean),by=.(day,run,stimNum,roi,soundID,winid)]
  out[,postStim:=winid>=0]
  return(out)
}
