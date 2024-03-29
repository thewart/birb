source("helperfuncs.R")
library(lme4)
library(data.table)
library(ggplot2)
path <- "/run/user/1001/gvfs/smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab/All_Staff/Seth/For Tom/"

brbid <- dir(path) |> str_subset("mat$",negate = T)
c <- c("Isolated","Female_distal_neg","Female_distal_pos","Female_proximal_neg","Female_proximal_pos")
nsamp <- 0
tsdat <- data.table()
ts0 <- vector()
shadup <- lmerControl(check.conv.singular = .makeCC(action = "ignore",tol = formals(isSingular)$tol))
for (b in brbid) {
  brbdat <- brbreader(file.path(path,b),cond=c,type=c("roi","stim"))
  
  for (i in 1:nrow(brbdat$meta)) {
    d <- brbdat$meta[i,day]
    r <- brbdat$meta[i,run]
    
    stim <- brbdat$stim[day==d & run==r]
    ldat <- brbdat$roi[day==d & run==r]
    ldat[,lum:=lm(log(lum) ~ 1 + t) |> resid(),by=roi]
    dummy <- stim[,soundFrames + shift(soundFrames)][-1]/2 |> round()
    stim <- rbind(stim,data.table(soundFrames=as.integer(dummy), stimNum=nrow(stim)+(1:length(dummy)),
                                  soundID=0, run=stim[1,run], day=stim[1,day]))
    
    redat <- binnify(ldat,stim,5,1)
    redat[,t := t/(15*60)]
    # redat[,postStim := as.numeric(postStim)]
    diffdat <- redat[,.(delta=diff(lum)),by=.(day,run,stimNum,soundID,roi)]
    # diffdat <- redat[,diff(lum),by=.(day,run,stimNum,soundID,roi)]
    
    # diffdat <- whiten_timepoints(ldat,stim,30,2,cov.est="shrinkage")
    
    brbdat
    lr <- logLik(lmer(delta ~ 1 + (1|roi) + (1|stimNum),data=diffdat,REML=F,control=shadup)) - 
      logLik(lmer(delta ~ 0 + (1|stimNum),data=diffdat,REML=F,control=shadup))
    tsdat <- rbind(tsdat,data.table(ts=lr,id=b,day=d,run=r,c=))
    
    if (nsamp>0) {
      tsi <- vector()
      ogd <- diffdat[,delta]
      for (j in 1:nsamp) {
        cat(j,"\r")
        diffdat[,delta:=permnull(ogd,stimNum)]
        # diffdat[,dlag:=shift(delta,fill=0),by=roi]
        tsi[j] <- logLik(lmer(delta ~ 1 + (1|roi) + (1|stimNum),data=diffdat,REML=F,control=shadup)) - 
          logLik(lmer(delta ~ 0 + (1|stimNum),data=diffdat,REML=F,control=shadup))
      }
      ts0 <- cbind(ts0,tsi)
      cat("\n")
    }
  }
}
