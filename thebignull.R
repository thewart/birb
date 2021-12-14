redat <- binnify(brbdat$roi[day==d],brbdat$stim[day==d],30,1)

nsamp <- 1e3
t0 <- matrix(ncol=nsamp,nrow=redat[,uniqueN(roi)])
for (i in 1:nsamp) {
  permdat <- redat[redat[,unique(postStim),by=.(t,stimNum)][,.(t,permStim=sample(V1)),by=stimNum],on=.(stimNum,t)]
  # permdat <- redat[,.(lum,permStim=sample(postStim,2)),by=.(run,stimNum,roi)]
  t0[,i] <- dcast(permdat,stimNum + roi ~ permStim, value.var="lum")[,t.test(log(`TRUE`)-log(`FALSE`))[[1]],by=roi]$V1
  # t0[i] <- dcast(permdat,stimNum + roi ~ permStim, value.var="lum")[,mean(`1`/`0`),by=roi][,median(V1-1)][[1]]
}
