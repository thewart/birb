dayreader <- function(path, cond="All", type=c("roi","stim")) {
  brbstr <- str_split_fixed(dir(path),"_",4)
  brbrun <- unique(brbstr[,3])
  filebase <- paste(dirname(path) |> basename(), basename(path),sep="_")
  roidat <- data.table()
  metadat <- data.table()
  stimdat <- data.table()
  
  for (r in brbrun) {
    brbpre <- file.path(path, paste(filebase,r,sep="_"))
    rmetadat <- fread(paste0(brbpre,"_info.csv"))
    if (!("All" %in% cond) & !(rmetadat$condition %in% cond)) next
    metadat <- rbind(metadat, rmetadat)
    
    if ("roi" %in% type) {
      rroi <- fread(paste0(brbpre,"_Fneuropil.csv")) |> t() |> as.data.table()
      rroi[,t:=1:.N]
      rroi <- melt(rroi, id.vars = "t", variable.name = "roi", value.name = "lum")
      # rroi[,lumz:=scale(lum),by=roi]
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

brbreader <- function(path, cond="All", type=c("roi","stim")) {
  brbday <- dir(file.path(path))
  out <- list(meta=data.table())
  for (l in type) out[[l]] <- data.table()
  
  for (d in brbday) {
    ddat <- dayreader(file.path(path,d), cond, type)
    if (length(ddat[[1]])>0) for (l in names(ddat)) {
      ddat[[l]]$day <- d
      out[[l]] <- rbind(out[[l]], ddat[[l]])
    }
  }
  
  return(out)
}
