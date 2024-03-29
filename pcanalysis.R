# c <- c("Baseline","Isolated","Female_distal_neg")
t0 <- 1
b <- mouseid[6]
brbdat <- brbreader(file.path(path,b))
ldat <- copy(brbdat$roi[t>t0])

# ldat[,lum:=lum-pil]
# ldat[,lum:=log(lum) - log(pil)]
ldat[,lum:=log(lum)-mean(lum),by=.(day,run,roi)]
ldat[,lum:=fffilt(lum,lo=4,padto = 2e4),by=.(day,run,roi)]
# ldat[,lum:=scale(lum),by=.(day,run,roi)]
# pcdat <- ldat[,(dcast(.SD, run + t ~ roi,value.var="lum") |> prcomp(center=T,scale=T,rank=5))$x |>
#                 cbind(t=min(t):max(t)) |> as.data.table() |> melt(id.vars="t"),by=.(day,run)]
pcdat <- ldat[,pcday(.SD,scale = F,rank=1),by=day]
pcdat <- ldat[,pcrun(.SD,scale = F,rank=1),by=.(day,run)]

pcdat <- pcdat[brbdat$beh[t>t0],on=.(day,run,t)]
pcdat <- pcdat[brbdat$meta[,.(run,day,condition)],on=.(day,run)]

pcdat[,cor(shoulder_energy,lum),by=.(day,run,component)] |> dcast(day ~ run)
pcdat[,.(abs(cor(shoulder_energy,lum)),mean(running_RPM>0)),by=.(day,run)][,qplot(y=V1,x=V2)]

cuts <- c(1,4,10,20,40,1000,2500,1e4)
lm(lum ~ fffilt(shoulder_energy,hi=4),data=pcdat[component=="PC1" & day=="220118" & run=="000"]) |> summary()
ggplot(pcdat[component=="PC1"],aes(x=t,y=lum)) + geom_line() +
  facet_grid(day ~ condition,scales="free")

melt(pcdat[day=="220128" & component=="PC1",.(run,t,lum,shoulder_energy)],id.vars = c("run","t")) |> 
  ggplot(aes(y=value,x=t)) + geom_line() + facet_grid(variable ~ run,scales = "free",space="free_x") + theme_classic()

pmove <- pcdat[day=="220128",density(shoulder_energy)]
mthresh <- pmove$x[match(2,diff(sign(diff(pmove$y))))+1]
