# mntdir <- "/run/user/1001/gvfs/"
mntdir <- "/run/user/1000/gvfs/"
path <- paste0(mntdir,"smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab/All_Staff/Seth/For Tom")
brbdir <- dir(path) |> str_subset("mat$",negate = T)

brbid <- brbdir[2]
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
    
    rroi <- fread(paste0(brbpre,"_Fneuropil.csv")) |> t() |> as.data.table()
    rroi[,t:=1:.N]
    rroi <- melt(rroi, id.vars = "t", variable.name = "roi")
    rroi[,z:=scale(value),by=roi]

    stimtim <- fread(paste0(brbpre,"_stim.csv"))
    rmetadata <- fread(paste0(brbpre,"_info.csv"))
    rmetadata$day <- brbday

    bs <- 5
    nb <- 1
    bindef <- data.table()
    for (t in 1:nrow(stimtim)) {
      st <- stimtim[t,soundFrames]
      trng <- st:(bs*nb+st-1)
      binid <- cut_number(trng,n = nb) |> as.numeric()
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

squishdat <- roidat[,.(t=mean(t)/(60*15),value=mean(value),z=mean(z),soundID=unique(soundID)),
                    by=.(roi,run,day,stimNum,binid)]
squishdat[,soundID:=as.factor(soundID)]
condat <- squishdat[metadat[condition=="USV_presentation"],on=.(run,day)]
condat[,postStim:=binid+1]
confit <- lmer(log(value) ~ 1 + day + t + postStim + (1|stimNum:day) + (1|stimNum:day:roi) + (1 + t|roi:day) +
                 (1 + postStim|roi), data=condat,REML=F)

baseformula <- "log(value) ~ 1 + day + t + postStim:soundID + (1|stimNum:day) + (1|stimNum:day:roi) + (1 + t|roi:day)"
soundReffs <- paste0("(0 + postStim:as.numeric(soundID==",condat[,unique(soundID)],")|roi:day)", collapse=" + ")

confit <- lmer(paste(baseformula,soundReffs,sep=" + "), data=condat,REML=F)

residuals(confit)[condat[,binid==0]] 

