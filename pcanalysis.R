# c <- c("Baseline","Isolated","Female_distal_neg")
t0 <- 1
b <- brbid[6]
brbdat <- brbreader(file.path(path,b))
ldat <- copy(brbdat$roi[t>t0])

# ldat[,lum:=lum-pil]
# ldat[,lum:=log(lum) - log(pil)]
ldat[,lum:=log(lum)]
ldat[,lum:=fffilt(lum-mean(lum),4,lopass=F,padto = 2e4),by=.(day,run,roi)]
# ldat[,lum:=scale(lum),by=.(day,run,roi)]
# pcdat <- ldat[,(dcast(.SD, run + t ~ roi,value.var="lum") |> prcomp(center=T,scale=T,rank=5))$x |>
#                 cbind(t=min(t):max(t)) |> as.data.table() |> melt(id.vars="t"),by=.(day,run)]
pcdat <- ldat[,pcday(.SD,scale = F),by=day]
pcdat <- ldat[,pcrun(.SD,scale = F),by=.(day,run)]

pcdat <- pcdat[brbdat$beh[t>t0],on=.(day,run,t)]
pcdat <- pcdat[brbdat$meta[,.(run,day,condition)],on=.(day,run)]

pcdat[,cor(shoulder_energy,lum),by=.(day,condition,component)] |> dcast(day + condition ~ component)
foo <- .Last.value$PC2
qplot(foo)

ggplot(pcdat[component=="PC1"],aes(x=t,y=lum)) + geom_line() +
  facet_grid(day ~ condition,scales="free")
