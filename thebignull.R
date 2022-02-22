b <- "TH135"
r <- "r007"
d <- "210615"
c <- "Baseline"
bdat <- brbreader(file.path(path,b),cond=c,type=c("roi","stim"))
# bdat$roi[,lumz := scale(lum),by=.(day,run,roi)]

stim <- bdat$stim[day==d & run==r]
# dummy <- stim[,soundFrames + shift(soundFrames)][-1]/2 |> round()
# stim <- rbind(stim,data.table(soundFrames=as.integer(dummy), stimNum=nrow(stim)+(1:length(dummy)), 
# soundID=0, run=stim[1,run], day=stim[1,day]))
ldat <- bdat$roi[day==d & run==r]
# udat <- whiten(ldat)
# lossybin(ldat,10)[,.(t=t[-1],delta=diff(log(lum))),by=.(day,run,roi)] |> dcast(t ~ roi, value.var="delta")

redat <- binnify(ldat,stim,30,1)
if (redat[,any(is.na(roi))]) redat <- redat[-which(is.na(roi))]
# uredat <- binnify(udat,stim[soundID==1],10,1)
# if (uredat[,any(is.na(roi))]) uredat <- uredat[-which(is.na(roi))]
# redat <- binnify(bdat$roi[day==d],bdat$stim[day==d],10,1)
# redat[,t := t/(15*60)]
redat[,postStim := as.numeric(postStim)]
uredat[,postStim := as.numeric(postStim)]

diffdat <- redat[,.(delta=diff(log(lum))),by=.(day,run,stimNum,soundID,roi)]
# udiffdat <- uredat[,.(delta=diff(lum)),by=.(day,run,stimNum,soundID,roi)]

# redat[,lum_norm:=qnorm((0:(.N+1))/(.N+1))[-c(1,.N+2)][rank(lum)],by=roi]
# diffdat <- redat[,diff(lum_norm),by=.(day,run,stimNum,soundID,roi)]

permnull(diffdat)

