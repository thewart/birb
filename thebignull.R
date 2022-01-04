d <- brbdat$meta[,unique(day)[3]]
redat <- binnify(brbdat$roi[day==d],brbdat$stim[day==d],30,1)

nsamp <- 1e3
t0 <- matrix(ncol=nsamp,nrow=redat[,uniqueN(roi)])
for (i in 1:nsamp) {
  permdat <- redat[redat[,unique(postStim),by=.(t,stimNum)][,.(t,permStim=sample(V1)),by=stimNum],on=.(stimNum,t)]
  # permdat <- redat[,.(lum,permStim=sample(postStim,2)),by=.(run,stimNum,roi)]
  t0[,i] <- dcast(permdat,stimNum + roi ~ permStim, value.var="lum")[,t.test(log(`TRUE`)-log(`FALSE`))[[1]],by=roi]$V1
  # t0[,i] <- dcast(permdat,stimNum + roi ~ permStim, value.var="lum")[,t.test(`TRUE`/`FALSE` - 1)[[1]],by=roi]$V1
}

#population mean test
apply(t0,2,function(x) t.test(x)[[1]]) |> abs() |> quantile(0.975)

#% significant test
apply(t0,2,function(x) mean(x>qt(0.975,df = nrow(t0)-1))) |> quantile(0.95)
