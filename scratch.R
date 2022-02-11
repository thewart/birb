brbid <- dir(path) |> str_subset("mat$",negate = T)
c <- "Baseline"
nsamp <- 1e3
lr <- vector()
ts <- vector()
for (b in brbid) {
  brbdat <- brbreader(file.path(path,b),cond=c,type=c("roi","stim"))
  
  for (i in 1:nrow(brbdat$meta)) {
    d <- brbdat$meta[i,day]
    r <- brbdat$meta[i,run]
    
    stim <- brbdat$stim[day==d & run==r]
    ldat <- brbdat$roi[day==d & run==r]
    # dummy <- stim[,soundFrames + shift(soundFrames)][-1]/2 |> round()
    # stim <- rbind(stim,data.table(soundFrames=as.integer(dummy), stimNum=nrow(stim)+(1:length(dummy)), 
    # soundID=0, run=stim[1,run], day=stim[1,day]))
    
    redat <- binnify(ldat,stim,10,1)
    redat[,t := t/(15*60)]
    # redat[,postStim := as.numeric(postStim)]
    diffdat <- redat[,.(delta=diff(log(lum))),by=.(day,run,stimNum,soundID,roi)]
    # diffdat <- redat[,diff(lum),by=.(day,run,stimNum,soundID,roi)]
    
    # diffdat <- whiten_timepoints(ldat,stim,30,2,cov.est="shrinkage")
    
    # lr <- rbind(lr,data.table(id = b, day = d, run = r,
    #                           whitened = logLik(lmer(z ~ 1 + (1|roi),data=diffdat,REML=F)) - 
    #                             logLik(lm(z ~ 0,data=diffdat,REML=F)),
    #                           rankone = logLik(lmer(delta ~ 1 + (1|roi) + (1|stimNum),data=diffdat,REML=F)) - 
    #                             logLik(lmer(delta ~ 0 + (1|stimNum),data=diffdat,REML=F)),
    #                           raw = logLik(lmer(delta ~ 1 + (1|roi),data=diffdat,REML=F)) - 
    #                             logLik(lm(delta ~ 0,data=diffdat,REML=F)))
    # )
    lr[i] <- logLik(lmer(delta ~ 1 + (1|roi) + (1|stimNum),data=diffdat,REML=F)) - 
      logLik(lmer(delta ~ 0 + (1|stimNum),data=diffdat,REML=F))
    
    tsi <- vector()
    ogd <- diffdat[,delta]
    for (j in 1:nsamp) {
      diffdat[,delta:=permnull(ogd,stimNum)]
      tsi[j] <- logLik(lmer(delta ~ 1 + (1|roi) + (1|stimNum),data=diffdat,REML=F)) - 
        logLik(lmer(delta ~ 0 + (1|stimNum),data=diffdat,REML=F))
    }
    ts <- cbind(ts,tsi)
  }
}

ggplot(melt(lr),aes(y=value,x=paste(id,day,run,sep = ":"),fill=variable)) + geom_col(position = position_dodge()) + 
  xlab(NULL) + ylab("Likelihood ratio") + scale_fill_discrete(NULL) + theme(legend.position = "bottom") +
  geom_hline(yintercept = qchisq(0.95,df = 2)/2)

##### dummy stims
stim <- bdat$stim[day==d]
dummy <- stim[,soundFrames + shift(soundFrames)][-1]/2 |> round()
stim <- rbind(stim,data.table(soundFrames=as.integer(dummy), stimNum=nrow(stim)+(1:length(dummy)), soundID=0, run=stim[1,run], day=stim[1,day]))
