ldat[,lum:=c(NA,diff(log(lum))),by=.(day,run,roi)]

c <- c("USV_presentation","USV_playback")
b <- brbid[1]
brbdat <- brbreader(file.path(path,b),cond = c)
ldat <- copy(brbdat$roi)[brbdat$beh,on=.(run,day,t)]

redat <- binnify(ldat,brbdat$stim,10,3)
# if not delta
redat[,lum:= log(lum) - log(lum[winid==-1]),by=.(day,run,roi,stimNum)]
redat <- redat[brbdat$meta[,.(run,day,condition)],on=.(run,day)]

tobs <- redat[winid>=0,.(lum=mean(lum),z=mean(lum)*sqrt(.N)/sd(lum)),by=.(day,condition,soundID,roi,winid)]
aveeff <- redat[winid>=0,mean(lum)*sqrt(.N)/sd(lum),by=.(day,condition,winid,roi)]
tobs[aveeff,on=.(day,condition,roi,winid)][V1>5 & winid==1,sd(lum)/mean(lum),by=.(day,condition,roi,winid)] |> 
  ggplot(aes(x=V1)) + geom_histogram() + facet_grid(condition ~ day) +
  ylab("# ROIs") + xlab("SD/mean luminance change") + ggtitle(paste(b,"variability of USV responses in active ROIs"))

# if not delta
# redat[,lum:= log(lum) - log(lum[winid==-1]),by=.(day,run,roi,stimNum)]
redat <- binnify(ldat,brbdat$stim,1,45)
redat[,lum:= log(lum) - log(lum[winid==-1]),by=.(day,run,roi,stimNum)]
redat <- redat[brbdat$meta[,.(run,day,condition)],on=.(run,day)]
redat[,mean(lum),by=.(day,condition,stimNum,soundID,winid)][
  ,mean(V1)*sqrt(.N)/sd(V1),by=.(day,condition,soundID,winid)] |> 
  ggplot(aes(y=V1,x=winid,color=factor(soundID))) + geom_line() + facet_grid(condition ~ day) 

redat <- binnify(ldat,brbdat$stim,10,3)
redat[,lum:= log(lum) - log(lum[winid==-1]),by=.(day,run,roi,stimNum)]
redat <- redat[brbdat$meta[,.(run,day,condition)],on=.(run,day)]
foo <- foo[,.(lum=mean(lum),rpm=mean(running_RPM),pupil=mean(pupil_pxlarea),whisker=mean(whisker_energy),
              nose=mean(nose_energy),jaw=mean(jaw_energy)),by=.(day,condition,soundID,winid,stimNum)]
foo <- foo[, `:=` (rpm=scale(rpm),pupil=scale(pupil),whisker=scale(whisker),nose=scale(nose),jaw=scale(jaw))]

lmer(lum ~ 1 + jaw + rpm + pupil + (1 + jaw + rpm + pupil|soundID),
     data=foo[day=="220118" & condition=="USV_playback"],REML=F) |> summary()

lmer(lum ~ 0 + day:condition + jaw + rpm + pupil + (1 + jaw + rpm + pupil|soundID:day:condition),data=foo,REML=F) |> summary()
lmer(lum ~ 1 + jaw + rpm + pupil + (1 + jaw + rpm + pupil|soundID:day:condition) + 
       (1 + jaw + rpm + pupil|day:condition),data=foo,REML=F) |> summary()



ldat[,lum:=log(lum)]
pcdat <- ldat[t>1,(dcast(.SD,t~roi,value.var="lum") |> prcomp(center=T,scale=T,rank=3))$x |>
                cbind(t=min(t):max(t)) |> as.data.table() |> melt(id.vars="t"),by=.(day,run)]
pcdat <- pcdat[brbdat$beh,on=.(day,run,t)]
setnames(pcdat,c("variable","value"),c("roi","lum"))
redat <- binnify(pcdat,brbdat$stim,30,1)
# if not delta
redat[,lum:= lum - lum[winid==-1],by=.(day,run,roi,stimNum)]
redat <- redat[brbdat$meta[,.(run,day,condition)],on=.(run,day)]
