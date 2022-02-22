# mntdir <- switch(Sys.info()['nodename'],
                 # `yawarakai-te`="/run/user/1001/gvfs",
                 # `seth-Z370P-D3`="/run/user/1000/gvfs")
# path <- file.path(mntdir,"smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab","All_Staff","Seth/For Tom")

brbid <- dir(path) |> str_subset("mat$",negate = T)
b <- brbid[2]
# brbday <- dir(file.path(path,b))
# d <- brbday[1]

c <- "Baseline"
brbdat <- brbreader(file.path(path,b),cond = c)
cmeta <- brbdat$meta
i <- 1
d <- cmeta[i,day]
r <- cmeta[i,run]

ggplot(brbdat$roi[day==d & run==r][roi %in% sample(unique(roi),6)][,.(t=t,lum=lum/mean(lum)),by=roi], aes(y=lum,x=t)) + 
  geom_line() + facet_wrap("roi") + ggtitle(paste(d,c))

redat <- binnify(brbdat$roi[day==d & run==r],brbdat$stim[day==d & run==r],15,1)
redat[,t := t/(15*60)]
redat[,postStim := as.numeric(postStim)]
# redat[,c("luml1","luml2") := .(c(NA,lum[-.N]),c(NA,NA,lum[-c(.N-1,.N)])),by=.(day,run,roi,stimNum)]

soundeffs <- paste0("postStim:as.numeric(soundID==",redat[,sort(unique(soundID))],")", collapse=" + ")
fixeffs <- paste0("log(lum) ~ 1 + t + ",soundeffs)
raneffs <- paste0("(1|t) + (1|stimNum:roi) + (1+t|roi) + (0+",soundeffs,"||roi)")

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

##### stan
model <- stan_model("caLciuMM_factor_intercept.stan")
standat <- list(Y = dcast(diffdat, roi ~ stimNum,value.var = "delta")[,-1,with=F] |> as.matrix(),
                N = diffdat[,uniqueN(stimNum)], M = diffdat[,uniqueN(roi)], K = 1)
fit <- sampling(model,standat,iter=400,chains=1,pars=c("mu_alpha","sigma_alpha","sigma_lambda","sigma","Lambda"))
