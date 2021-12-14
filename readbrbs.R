mntdir <- switch(Sys.info()['nodename'],
                 `yawarakai-te`="/run/user/1001/gvfs/",
                 `seth-Z370P-D3`="/run/user/1000/gvfs/")
path <- file.path(mntdir,"smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab","All_Staff","Seth/For Tom")

brbid <- dir(path) |> str_subset("mat$",negate = T)
b <- brbid[1]
# brbday <- dir(file.path(path,b))

brbdat <- brbreader(file.path(path,b),cond = "USV_presentation")
d <- brbdat$meta[,unique(day)[1]]
# daydat <- dayreader(file.path(path,b,d), cond = "Baseline")

redat <- binnify(brbdat$roi[day==d],brbdat$stim[day==d],30,1)

confit <- lmer(log(value) ~ 1 + t + postStim + (1|stimNum) + (1|stimNum:roi) + (1 + t|roi) +
                 (0 + postStim|roi), data=condat,REML=F)

baseformula <- "log(pil) ~ 1 + t + postStim:as.factor(soundID) + (1|stimNum) + (1|stimNum:roi) + (1 + t|roi)"
soundReffs <- paste0("(0 + postStim:as.numeric(soundID==",condat[,unique(soundID)],")|roi)", collapse=" + ")

confit <- lmer(paste(baseformula,soundReffs,sep=" + "), data=condat,REML=F)
bigfit <- lmer(paste(baseformula,"(1|t)",soundReffs,sep=" + "), data=condat,REML=F)

residuals(confit)[condat[,binid==0]] 



