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
