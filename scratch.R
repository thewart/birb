# mntdir <- switch(Sys.info()['nodename'],
# `yawarakai-te`="/run/user/1001/gvfs",
# `seth-Z370P-D3`="/run/user/1000/gvfs")
# path <- file.path(mntdir,"smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab","All_Staff","Seth/For Tom")

brbid <- dir(path) |> str_subset("mat$",negate = T)
c <- c("Female_distal_neg","Female_distal_pos","Female_proximal_neg","Female_proximal_pos")
b <- brbid[6]
brbdat <- brbreader(file.path(path,b),cond = c)


d <- brbdat$meta[,unique(day)][1]
stim <- brbdat$stim[day==d]
ldat <- brbdat$roi[day==d]
ldat[,lum:=lm(log(lum) ~ 1 + t) |> resid(),by=.(run,roi)]
# dummy <- stim[,soundFrames + shift(soundFrames)][-1]/2 |> round()
# stim <- rbind(stim,data.table(soundFrames=as.integer(dummy), stimNum=nrow(stim)+(1:length(dummy)),
#                               soundID=0, run=stim[1,run], day=stim[1,day]))
# setkey(stim,"soundFrames")
# stim[,stimNum:=1:.N]

redat <- binnify(ldat,stim,15,1)
redat[,t := t/(15*60)]
# redat[,postStim := as.numeric(postStim)]
diffdat <- redat[,.(delta=diff(lum)),by=.(run,stimNum,soundID,roi)]

ggplot(brbdat$roi[day==d & run==r][roi %in% sample(unique(roi),6)][,.(t=t,lum=lum/mean(lum)),by=roi], aes(y=lum,x=t)) + 
  geom_line() + facet_wrap("roi") + ggtitle(paste(d,c))

soundeffs <- paste0("as.numeric(soundID==",redat[,sort(unique(soundID))],")")
fixeffs <- paste0("log(lum) ~ 1 + ",soundeffs)
raneffs <- paste0("(1|stimNum:run) + (0+",soundeffs,"||roi)")

bigfit <- lmer(paste(fixeffs,raneffs,sep=" + "), data=redat,REML=F)
lilfit <- lmer(log(lum) ~ 1 + t + postStim + (1|t) + (1|stimNum:roi) + (1+t|roi) + 
                 (0+postStim|soundID) + (0+postStim|roi) + (0+postStim|soundID:roi), data=redat,REML=F)
nullfit <- lmer(log(lum) ~ 1 + t + (1|t) + (1|stimNum:roi) + (1+t|roi), data=redat,REML=F)

anova(bigfit,lilfit,nullfit)


##### single-roi time series analysis w/ behavior 
foo <- daydat$roi[daydat$beh,on=.(run,t)]
foo[,sec:=cut_width(t,15,7.5) |> as.numeric()]
foo <- foo[,lapply(.SD,mean),by=.(run,roi,sec)]
lm(log(lum) ~ running_RPM + shift(log(lum),1) + shift(running_RPM,1) + shift(log(lum),2) + shift(running_RMP,2),data=foo[run=="r001" & roi=="V1"])

ggplot(melt(lr),aes(y=value,x=paste(id,day,run,sep = ":"),fill=variable)) + geom_col(position = position_dodge()) + 
  xlab(NULL) + ylab("Likelihood ratio") + scale_fill_discrete(NULL) + theme(legend.position = "bottom") +
  geom_hline(yintercept = qchisq(0.95,df = 2)/2)

##### dummy stims
stim <- bdat$stim[day==d]
dummy <- stim[,soundFrames + shift(soundFrames)][-1]/2 |> round()
stim <- rbind(stim,data.table(soundFrames=as.integer(dummy), stimNum=nrow(stim)+(1:length(dummy)), soundID=0, run=stim[1,run], day=stim[1,day]))
setkey(stim,"soundFrames")
stim[,stimNum:=1:.N]

##### stan
model <- stan_model("caLciuMM_factor_intercept.stan")
standat <- list(Y = dcast(diffdat, roi ~ stimNum,value.var = "delta")[,-1,with=F] |> as.matrix(),
                N = diffdat[,uniqueN(stimNum)], M = diffdat[,uniqueN(roi)], K = 1)
fit <- sampling(model,standat,iter=400,chains=1,pars=c("mu_alpha","sigma_alpha","sigma_lambda","sigma","Lambda"))


##### psth
ldat <- brbdat$roi
ldat[,lum:=lm(log(lum) ~ 1 + t) |> resid(),by=.(roi,run)]
redat <- binnify(ldat,brbdat$stim,1,15)
diffdat <- redat[,.(kHz=soundFrequency[-1]/1e3,winid=winid[-1],delta=diff(lum)),
                 by=.(day,run,stimNum,roi)]

diffdat[winid>-6,mean(delta),by=.(day,run,kHz,winid)][brbdat$meta[,.(day,run,condition)],on=.(day,run)] |> 
  ggplot(aes(x=winid,y=V1,color=ordered(kHz))) + geom_line() + facet_grid(day ~ condition) +
  geom_hline(yintercept = 0,color="gray") + geom_vline(xintercept = 0,color="gray") + 
  ylab("delta") + xlab("timepoint") + ggtitle(b) + theme_classic()

shadup <- lmerControl(check.conv.singular = .makeCC(action = "ignore",tol = formals(isSingular)$tol))

lldat <- data.table()
manyfits <- list()
dses <- diffdat[,.(day,kHz)] |> unique()
for (d in dses[,unique(day)]) {
  # manyfits[[d]] <- list()
  ses <- dses[day==d,kHz] |> sort()
  
  for (s in ses) {
    # manyfits[[d]][[paste0(s,"kHz")]] <- list()
    
    for (t in -5:14) {
      
      cat(paste(d,s,t,sep=":"),"\r")
      
      dstdat <- diffdat[day==d & kHz==s & winid==t,.(run,roi,stimNum,delta)]
      
      fit1 <- lmer(delta ~ 1 + (1|roi) + (1|stimNum:run),data=dstdat,REML=F,control=shadup)
      fit0 <- lmer(delta ~ 0 + (1|stimNum:run),data=dstdat,REML=F,control=shadup)
      
      if (dstdat[,uniqueN(run)] > 1) {
        dummified <- model.matrix(~ 0 + run,dstdat)
        dstdat <- cbind(dstdat,dummified)
        runs <- paste(colnames(dummified),collapse =" + ")
        # h2 <- paste("delta ~ 0", runs, "(1 | roi)", paste0("(0 + ",runs," || roi)"),
        #             "(1 | stimNum:run)",sep = " + ")
        # fit2 <- lmer(as.formula(h2),data=dstdat,REML=F,control=shadup)
        fit2 <- lmer(delta ~ 1 + run + (1|roi) + (1|roi:run) + (1|stimNum:run),data=dstdat,REML=F,control=shadup)
        manyfits[[d]][[paste0(s,"kHz")]][[paste0("t",t)]] <- fit2
        ll2 <- logLik(fit2) |> as.numeric()
        
      } else {
        
        ll2 <- NA_real_
      }
      
      lldat <- rbind(lldat,data.table(ll0=logLik(fit0),ll1=logLik(fit1),ll2=as.numeric(ll2),day=d,kHz=s,t=t))
    }
  }
}

