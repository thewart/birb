fffilt <- function(og, lo=1, hi=padto/2, padto=2e4) {
  npad <- padto - length(og)
  og <- c(og, rep(mean(og),npad))
  juh <- fft(og-mean(og))
  
  filt <- rep(0,padto/2)
  filt[lo:hi] <- 1
  filt <- c(0,filt,rev(filt)[-1])

  return( (Re(fft(filt*juh,inverse = T))/length(juh))[1:(length(juh)-npad)] + mean(og))
}



# brbdat$roi <- brbdat$roi[,bl:=lopass(lum,10),by=.(day,roi)]
# condat[,postStim:=as.numeric(winid>=0)]
# 
# lmer(log(lum) ~ 1 + offset(log(bl)) + t + postStim + (1|t) + (1|stimNum) + (1|stimNum:roi) + (1 + t|roi) + (0 + postStim|roi), data=condat,REML=F) |> summary()
# 
# ggplot(brbdat$roi[day==unique(day)[1] & roi=="V3"]) + geom_line(aes(x=t,y=lum)) + geom_line(aes(x=t,y=bl),color="red")
