lopass <- function(og, lofi) {
  juh <- fft(og-mean(og))
  filt <- c(rep(1,lofi),rep(0,length(juh)-lofi*2),rep(1,lofi))
  
  return(Re(fft(filt*juh,inverse = T))/length(juh) + mean(og))
}

filtdat <- lildat[lildat[,.(t,bl=lopass(value,300)),by=roi],on=.(roi,t)]
condat <- filtdat[stimdat[,slicetime(soundFrames,20,1),by=.(day,run,stimNum,soundID)][
  metadat[condition=="Baseline",],on=.(day,run)][day=="210511"],
  on=.(run,t)][,.(t=mean(t)/(60*15),lum=mean(value),bl=mean(bl)),
                   by=.(day,run,stimNum,roi,soundID,winid)]
condat[,postStim:=as.numeric(winid>=0)]

lmer(log(lum) ~ 1 + offset(log(bl)) + t + postStim + (1|t) + (1|stimNum) + (1|stimNum:roi) + (1 + t|roi) + (0 + postStim|roi), data=condat,REML=F) |> summary()

ggplot(filtdat[roi=="V138"]) + geom_line(aes(x=t,y=value)) + geom_line(aes(x=t,y=bl),color="red")
