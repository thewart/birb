# mntdir <- "/run/user/1001/gvfs/"
mntdir <- "/run/user/1000/gvfs/"
path <- paste0(mntdir,"smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab/All_Staff/Seth/For Tom")
brbdir <- dir(path) |> str_subset("mat$",negate = T)

brbid <- brbdir[1]
brblife <- dir(paste(path,brbid,sep="/"))

roidat <- data.table()
metadat <- data.table()

for (d in 1:length(brblife)) {
  brbday <- brblife[d]
  brbdat <- dir(paste(path,brbid,brbday,sep="/"))
  brbstr <- str_split_fixed(brbdat,"_",4)
  brbrun <- unique(brbstr[,3])
  
  for (r in 1:length(brbrun)) {
    brbpre <- paste0(paste(path,brbid,brbday,"/",sep="/"), paste(brbid,brbday,brbrun[r],sep="_"))
    
    rroi <- fread(paste0(brbpre,"_ROIs.csv")) |> t() |> as.data.table()
    rroi[,t:=1:.N]
    rroi <- melt(rroi, id.vars = "t", variable.name = "roi")
    rroi[,value := scale(value),by=roi]

    stimtim <- fread(paste0(brbpre,"_stim.csv"))
    rmetadata <- fread(paste0(brbpre,"_info.csv"))
    rmetadata$day <- brbday

    bs <- 5
    nb <- 1
    bindef <- data.table()
    for (t in 1:nrow(stimtim)) {
      st <- stimtim[t,soundFrames]
      trng <- st:(bs*nb+st-1)
      binid <- cut_number(trng,n = nb) %>% as.numeric()
      binid <- c(sort(-binid),binid-1)
      trng <- c((st-bs*nb):(st-1),trng)
      bindef <- rbind(bindef, data.table(t=trng,binid=binid,
                                         soundID=rep(stimtim[t,soundID],length(binid)),
                                         stimNum=rep(t,length(binid))))
    }
    rroi <- rroi[bindef,on="t"]
    rroi$run <- rmetadata$run
    rroi$day <- rmetadata$day
    
    roidat <- rbind(roidat, rroi)
    metadat <- rbind(metadat, rmetadata)
    
  }
}

