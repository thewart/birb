path <- "/run/user/1001/gvfs/smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab/All_Staff/Seth/For Tom/"
source("helperfuncs.R")
mouseid <- dir(path)

thesevars <- c("von_onset","shoulder_delta","rpm_delta") # variables in the full model
nothesevars <- c("shoulder_delta","rpm_delta") # variables in the null model
fitdat <- data.table()
sigdat <- data.table()
for (thismouse in mouseid) {
  # thismouse <- mouseid[6]
  behdat <- brbreader(file.path(path,thismouse),type="beh")
  
  vocsess <- behdat$beh[,sum(vocal_present),by=.(day,run)][V1>0]
  vocsess <- vocsess[day %in% vocsess[,sum(V1),by=day][V1>900,day]]
  vocsess <- behdat$meta[vocsess,on=.(day,run)][,.(day,run,condition)]
  vocbeh <- behdat$beh[vocsess,on=.(day,run)]

  for (thisday in vocsess[,unique(day)]) {
    # thisday <- vocsess[,unique(day)[j]]
    theseruns <- vocsess[day==thisday,run]
    mousedat <- dayreader(file.path(path,thismouse,thisday),run=theseruns)
    
    # vocbouts <- vocbeh[day==thisday,get_bouts(vocal_present),by=.(condition,run)]
    # vocwindows <- vocbouts[,c(dur_bout=dur_bout,slicetime(start_bout,1,15)),by=.(condition,run,num_bout)]
    # vocdat <- ldat[vocwindows,on=.(run,t)]
    
    # ggplot(vocdat[,mean(lum),by=.(winid,day)],aes(x=winid,y=V1,color=day)) + geom_line()
    # ggplot(vocdat[,mean(lum),by=.(winid,roi)],aes(x=winid,y=V1,group=roi)) + geom_line(alpha=0.1)
    
    delta_dat <- vocbeh[,.(shoulder_delta=c(NA,diff(shoulder_deltapix)),
                           rpm_delta=c(NA,diff(running_RPM)),
                           vrate_delta=c(NA,diff(vocal_rate)),
                           von_delta=c(NA,diff(vocal_present)),
                           bp_delta=c(NA,diff(bandpower_USV)),
                           fdist10_delta=c(NA,diff(Female_distance<10)),
                           fdist5_delta=c(NA,diff(Female_distance<5))),
                        by=.(day,run)]
    delta_dat[,`:=` (von_onset=as.numeric(von_delta==1),
                     von_offset=as.numeric(von_delta==-1),
                     fdist10_onset=as.numeric(fdist10_delta==1),
                     fdist10_offset=as.numeric(fdist10_delta==-1),
                     fdist10_onset=as.numeric(fdist10_delta==1),
                     fdist10_offset=as.numeric(fdist10_delta==-1)),
              by=.(day,run)]
                     
                     
    delta_dat <- cbind(vocbeh[,.(t)],delta_dat[,c(
      shift(shoulder_delta,-5:15,give.names=T),
      shift(vrate_delta,-5:15,give.names = T),
      shift(rpm_delta,-5:15,give.names=T),
      shift(von_onset,-5:15,give.names=T),
      shift(von_offset,-5:15,give.names=T),
      shift(bp_delta,-5:15,give.names=T)),
      # shift(von_delta,-5:15,give.names=T)),
      by=.(day,run)])
    
    ldat <- mousedat$roi
    # ldat[,loglum:=fffilt(lum, lo=4, padto = 2.5e4),by=.(run,roi)]
    ldat[,loglum := log(lum) - mean(log(lum)), by= .(run,roi)]
    ldat[,loglum_delta := c(NA,diff(loglum)),by=.(run,roi)]
    # ldat[,logpil_delta := c(NA,diff(log(pil))),by=.(run,roi)]
    ldat <- cbind(ldat,ldat[,shift(loglum_delta,1:15,give.names = T),by=.(run,roi)][,-c(1,2),with=F])
    # ldat <- cbind(ldat,ldat[,shift(loglum,1:15,give.names = T),by=.(run,roi)][,-c(1,2),with=F])
    
    for (i in 1:ldat[,uniqueN(roi)]) {
      cat(thismouse, " ", thisday, " ", i,"/",ldat[,uniqueN(roi)], "\r")
      
      covariates <- delta_dat[day==thisday,c("run","t",sapply(thesevars,str_subset,string=names(delta_dat)) |> as.vector()),with=F]
      fitme <- ldat[roi==paste0("V",i)][covariates,on=.(run,t)][,-c(1:6),with=F]
      nofitme <- fitme[!is.na(rowSums(fitme))]
      fullmodel <- lm(loglum_delta ~ ., data=nofitme)
      nullmodel <- lm(loglum_delta ~ ., data=nofitme[,sapply(c("loglum_delta",nothesevars),str_subset,string=names(nofitme)) |> unlist(),with=F])
      cov_info <- str_split_fixed(names(fullmodel$coefficients[-1]),"_",4)

      fitdat <- rbind(fitdat,
                      data.table(beta = fullmodel$coefficients[-1],
                                 tstat = fullmodel$coefficients[-1]/sqrt(diag(vcov(fullmodel))[-1]),
                                 covariate = apply(cov_info[,1:2],1,paste0,collapse="_"), 
                                 lag = as.numeric(cov_info[,4]) * ifelse(cov_info[,3]=="lead",-1,1), 
                                 roi = paste0("V",i),day = thisday, mouse = thismouse))
      sigdat <- rbind(sigdat,
                      data.table(p = anova(nullmodel,fullmodel)$`Pr(>F)`[2],
                                 roi = paste0("V",i), day = thisday, mouse = thismouse))
    }
    cat("\n")
  }
}

sigunits <- fitdat[covariate == "vocal_onset"][sigdat[p<.01],on=.(day,roi)]
wideunits <- dcast(sigunits,formula = roi + day ~ covariate + lag, value.var="tstat")[,-"roi",with=F]

# prout <-prcomp(X,scale=T)
# as.data.table(melt(prout$rotation[,1:3]))[,x:=1:.N,by=Var2][,c("cov","lag") := data.table(str_split_fixed(Var1,"_",2))] |> 
#   ggplot(aes(y=value,x=as.numeric(lag),color=Var2)) + geom_line() + facet_wrap(vars(cov),nrow=1,scales="free_x") + 
#   geom_hline(yintercept = 0) + theme_cowplot()

sigdat[,pfdr:=p.adjust(p,method="fdr"),by=.(mouse,day)]
ggplot(sigdat,aes(x=p)) + geom_histogram() + facet_wrap(vars(mouse),ncol=5)
fitdat[covariate=="von_onset"][sigdat[pfdr<0.05],on=.(day,mouse,roi)] |> 
  ggplot(aes(y=beta,x=lag,group=roi,color=mouse)) + geom_line(alpha=0.25) + facet_wrap(vars(paste(mouse,day,sep=":")),scales="free") + theme_minimal()
fitdat[covariate=="von_onset"][sigdat[pfdr<0.05],on=.(day,mouse,roi)][day %in% sigdat[pfdr<0.05][,.N,by=.(day,mouse)][N>5,day]][,mean(tstat),by=.(mouse,day,lag)] |> 
  ggplot(aes(y=V1,x=lag,color=mouse)) + geom_line() + facet_wrap(vars(paste(mouse,day,sep=":")),scales="free") + theme_minimal()

fitdat[covariate=="shoulder_delta"] |> ggplot(aes(y=beta,x=lag,group=roi,color=mouse)) + 
  geom_line(alpha=0.1) + facet_wrap(vars(paste(mouse,day,sep=":")),scales="free") + theme_minimal()
