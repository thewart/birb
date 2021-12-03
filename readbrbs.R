mntdir <- switch(Sys.info()['nodename'],
                 `yawarakai-te`="/run/user/1001/gvfs/",
                 `seth-Z370P-D3`="/run/user/1000/gvfs/")
path <- file.path(mntdir,"smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab","All_Staff","Seth/For Tom")
brbdir <- dir(path) |> str_subset("mat$",negate = T)

brbid <- brbdir[2]
brblife <- dir(file.path(path,brbid))
cond <- c("Baseline")

roidat <- data.table()
metadat <- data.table()
stimdat <- data.table()


confit <- lmer(log(value) ~ 1 + day + t + postStim + (1|stimNum:day) + (1|stimNum:day:roi) + (1 + t|roi:day) +
                 (0 + postStim|roi:day), data=condat,REML=F)

baseformula <- "log(lum) ~ 1 + day + t + postStim:soundID + (1|stimNum:day) + (1|stimNum:day:roi) + (1 + t|roi:day)"
soundReffs <- paste0("(0 + postStim:as.numeric(soundID==",condat[,unique(soundID)],")|roi:day)", collapse=" + ")

confit <- lmer(paste(baseformula,soundReffs,sep=" + "), data=condat,REML=F)

residuals(confit)[condat[,binid==0]] 

slicetime <- function(t,wsize,wnum) {
  trng <- t:(wsize*wnum+t-1)
  winid <- cut_number(trng,n = wnum) |> as.numeric()
  winid <- c(sort(-winid),winid-1)
  trng <- c((t-wsize*wnum):(t-1),trng)
  return(data.table(t=trng,winid=winid))
}


condat <- roidat[stimdat[,slicetime(soundFrames,20,1),by=.(day,run,stimNum,soundID)][
  metadat[condition=="Baseline",],on=.(day,run)],
  on=.(day,run,t)][,.(t=mean(t)/(60*15),lum=mean(value),lumz=mean(z)),
                   by=.(day,run,stimNum,roi,soundID,winid)]
condat[,postStim:=as.numeric(winid>=0)]
