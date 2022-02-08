brbid <- dir(path) |> str_subset("mat$",negate = T)
c <- "Baseline"
lr <- data.table()
for (b in brbid) {
  brbdat <- brbreader(file.path(path,b),cond=c,type=c("roi","stim"))
  
  for (j in 1:nrow(brbdat$meta)) {
    d <- brbdat$meta[j,day]
    r <- brbdat$meta[j,run]
    
    stim <- brbdat$stim[day==d & run==r]
    ldat <- brbdat$roi[day==d & run==r]
    # udat <- whiten(ldat[t != 1])
    # dummy <- stim[,soundFrames + shift(soundFrames)][-1]/2 |> round()
    # stim <- rbind(stim,data.table(soundFrames=as.integer(dummy), stimNum=nrow(stim)+(1:length(dummy)), 
    # soundID=0, run=stim[1,run], day=stim[1,day]))
    
    # redat <- binnify(udat,stim,30,1)
    # redat[,t := t/(15*60)]
    # redat[,postStim := as.numeric(postStim)]
    # diffdat <- redat[,diff(log(lum)),by=.(day,run,stimNum,soundID,roi)]
    # diffdat <- redat[,diff(lum),by=.(day,run,stimNum,soundID,roi)]
    
    diffdat <- whiten_timepoints(ldat,stim,30,2,cov.est="shrinkage")
    
    lr <- rbind(lr,data.table(id = b, day = d, run = r,
                              whitened = logLik(lmer(z ~ 1 + (1|roi),data=diffdat,REML=F)) - 
                                logLik(lm(z ~ 0,data=diffdat,REML=F)),
                              rankone = logLik(lmer(delta ~ 1 + (1|roi) + (1|stimNum),data=diffdat,REML=F)) - 
                                logLik(lmer(delta ~ 0 + (1|stimNum),data=diffdat,REML=F)),
                              raw = logLik(lmer(delta ~ 1 + (1|roi),data=diffdat,REML=F)) - 
                                logLik(lm(delta ~ 0,data=diffdat,REML=F)))
    )
  }
}

ggplot(melt(lr),aes(y=value,x=paste(id,day,run,sep = ":"),fill=variable)) + geom_col(position = position_dodge()) + 
  xlab(NULL) + ylab("Likelihood ratio") + scale_fill_discrete(NULL) + theme(legend.position = "bottom") +
  geom_hline(yintercept = qchisq(0.95,df = 2)/2)

##### dummy stims
stim <- bdat$stim[day==d]
dummy <- stim[,soundFrames + shift(soundFrames)][-1]/2 |> round()
stim <- rbind(stim,data.table(soundFrames=as.integer(dummy), stimNum=nrow(stim)+(1:length(dummy)), soundID=0, run=stim[1,run], day=stim[1,day]))
