path <- "/run/user/1001/gvfs/smb-share:server=duhsnas-pri.dhe.duke.edu,share=dusom_mooneylab/All_Staff/Seth/For Tom"
brbdir <- dir(path) |> str_subset("mat$",negate = T)

brbid <- brbdir[1]
brblife <- dir(paste(path,brbid,sep="/"))
brbday <- brblife[1]
brbdat <- dir(paste(path,brbid,brbday,sep="/"))
brbstr <- str_split_fixed(brbdat,"_",4)
brbrun <- unique(brbstr[,3])

brbpre <- paste0(paste(path,brbid,brbday,"/",sep="/"), paste(brbid,brbday,brbrun[1],sep="_"))

roi <- fread(paste0(brbpre,"_ROIs.csv")) |> t() |> as.data.table()
roi[,t:=1:.N]
roi <- melt(roi, id.vars = "t", variable.name = "roi")

stimtim <- fread(paste0(brbpre,"_stim.csv"))

bs <- 2
nb <- 5

bindef <- data.table()
for (i in 1:nrow(stimtim)) {
  st <- stimtim[i,soundFrames]
  trng <- st:(bs*nb+st-1)
  binid <- cut_number(trng,n = nb) %>% as.numeric()
  binid <- c(sort(-binid),binid-1)
  bindef <- rbind(bindef, data.table(t=trng,binid=binid))
}
