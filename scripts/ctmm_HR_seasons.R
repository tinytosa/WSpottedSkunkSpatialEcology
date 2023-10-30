##########################
#explore seasonality of home ranges

##########################

require(lubridate)

#######################
#load data
hr_ca_skunks <- read.table(file="data/HR_ctmm/hr_ca_estimates.txt", sep=",", header=T)

f <- c("SG-005", "SG-007", "SG-008", "SG-009", "SG-015", "SG-016", "SG-017","SG-018","SG-020") #females
nad83z10 <- "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

load("data/HR_ctmm/DATA_all_skunks_error_2022-11-22.rda") #DATA file for all skunks with error ellipses (vhf + gps)
load("data/HR_ctmm/FITS_all_skunks_error_2022-11-16.rda") #FITS file

skunk.info <- read.table("data/AllSkunkCapts.txt", sep=",", header=T)
skunk.info.unique <- ddply(skunk.info, .(Animal.ID, Sex), summarize, 
                           mass=mean(Mass..g., na.rm=T), mass.min=min(Mass..g., na.rm=T), mass.max=max(Mass..g., na.rm=T),
                           length.total=mean(Length.Total..cm., na.rm=T))

skunk.info.unique %>% group_by(Sex) %>% dplyr::summarize(mean = mean(mass, rm.na=T), n=length(mass), sd=sd(mass), se = sd(mass)/sqrt(length(mass)))
skunk.info.unique %>% group_by(Sex) %>% dplyr::summarize(mean = mean(length.total, rm.na=T), n=length(length.total), sd=sd(length.total), se = sd(length.total)/sqrt(length(length.total)))

skunk.info.unique <-  rbind(skunk.info.unique, data.frame(Animal.ID = "SG-001_2", Sex="M",
                                                          mass=595, mass.min=550, mass.max=660,
                                                          length.total=41.1))

#######################
# m1: `SG-002`, `SG-006`, `SG-011`, `SG-013`, `SG-021`, `SG-023`, `SG-025`
# m2: `SG-001`, `SG-001_2`, `SG-003`, `SG-004`, `SG-010`, `SG-012`, `SG-014`, `SG-019`, `SG-022`
# f:  `SG-005`, `SG-007`, `SG-008`, `SG-009`,`SG-015`,`SG-016`,`SG-017`,`SG-018`,`SG-020`
#######################

#seasons
#summer:  June - September
#fall:    October - January
#spring:  February - May

data <- NULL

for(i in 1:length(DATA))
{
  print(names(DATA)[i])
  data <- bind_rows(data, data.frame(animalid=names(DATA)[i], DATA[[i]]))
}

data <- merge(data, hr_ca_skunks[hr_ca_skunks$size=="homerange",c("animal","group")], by.x="animalid", by.y="animal", all.x=T)

p <- ggplot(data) + geom_point(aes(x=timestamp, y=animalid, col=group)) + xlab("date") + ylab("") + theme_bw(base_size=20)
p
# ggsave(p, filename="Figures/location_timestamps.tiff", height=6.5, width=10, units="in", dpi=330, compression="lzw")


data$month <- month(data$timestamp)
data$season <- "fall"
data[data$month > 5 & data$month < 10,]$season <- "summer"
data[data$month > 1 & data$month < 6,]$season <- "spring"

m <- data.frame(table(data[,c("animalid","month")]))
m <- merge(m, hr_ca_skunks[,c("animal","group")], by.x="animalid", by.y="animal")

table(data[,c("animalid","month")])
table(data[,c("animalid","season")])

num.locs <- data.frame(table(data[,c("animalid","season")]))
m[m$Freq == 0,]$group <- NA 
ggplot(m) + geom_tile(aes(x=month, y=animalid, fill=group)) + theme_bw(base_size=20)


###########################################################

sk.locs <- read.table("data/HR_ctmm/as_telemetry_obj_vhf_2022-11-17.txt", sep=",", header=T)
sk.locs$y <- year(sk.locs$timestamp)

sk.spring <- sk.locs[month(sk.locs$timestamp) < 6 & month(sk.locs$timestamp) > 1,]
sk.spring <- sk.spring[!sk.spring$individual.local.identifier %in% c("SG-002","SG-014","SG-015","SG-017","SG-019","SG-023"),]
table(sk.spring[,c("individual.local.identifier","y")])

sk.summer <- sk.locs[month(sk.locs$timestamp) < 10 & month(sk.locs$timestamp) > 5,]
sk.summer <- sk.summer[!sk.summer$individual.local.identifier %in% c("SG-012","SG-014","SG-016","SG-017","SG-023","SG-025"),]
table(sk.summer[,c("individual.local.identifier","y")])

sk.fall <- sk.locs[month(sk.locs$timestamp) > 9 | month(sk.locs$timestamp) == 1,]
table(sk.fall[,c("individual.local.identifier","y")])

###########################################################
#create DATA file for each season and save

# sk.locs <- sk.spring
# sk.locs <- sk.summer
# sk.locs <- sk.fall
# DATA <- as.telemetry(sk.locs, projection=CRS(nad83z10))
# 
# uere(DATA) <- 1 # you need this in newer versions of the package!
# 
# #read in error ellipse information for each skunk
# for(i in 1:length(DATA))
# {
#   animal <- names(DATA)[[i]]
#   print(animal)
#   
#   DATA[[i]]$COV.x.x <- sk.locs[sk.locs$individual.local.identifier == animal,]$COV.x.x
#   DATA[[i]]$COV.y.y <- sk.locs[sk.locs$individual.local.identifier == animal,]$COV.y.y
#   DATA[[i]]$COV.x.y <- sk.locs[sk.locs$individual.local.identifier == animal,]$COV.x.y
#   DATA[[i]]$VAR.xy <- (DATA[[i]]$COV.x.x + DATA[[i]]$COV.y.y)/2
# }
# 
# #plot vhf locations with error ellipses
# COL <- rainbow(length(DATA))
# plot(DATA, col=COL, lwd=2)
# 
# # SPRING <- DATA
# names(SPRING) <- paste(names(SPRING), ".spring", sep="")
# save(SPRING, file="data/HR_ctmm/DATA-SPRING_all_skunks_error_2023-01-03.rda")
# # SUMMER <- DATA
# names(SUMMER) <- paste(names(SUMMER), ".summer", sep="")
# save(SUMMER, file="data/HR_ctmm/DATA-SUMMER_all_skunks_error_2023-01-03.rda")
# # FALL <- DATA
# names(FALL) <- paste(names(FALL), ".fall", sep="")
# save(FALL, file="data/HR_ctmm/DATA-FALL_all_skunks_error_2023-01-03.rda")

###########################################################

load("data/HR_ctmm/DATA-SPRING_all_skunks_error_2023-01-03.rda")
load("data/HR_ctmm/DATA-SUMMER_all_skunks_error_2023-01-03.rda")
load("data/HR_ctmm/DATA-FALL_all_skunks_error_2023-01-03.rda")

DATA <- c(SPRING, SUMMER, FALL) #combine spring, summer, and fall lists

##########################################################
#skip this and load
FITS <- list()
for(i in 1:length(DATA))
{
  animal <- names(DATA)[[i]]
  print(animal)
  COL <- color(DATA[[i]], by='time')
  GUESS <- ctmm.guess(DATA[[i]], interactive=F, CTMM = ctmm(error=T)) #CTMM-ctmm(error=T) include this 
  FITS[[i]] <- ctmm.select(DATA[[i]], GUESS, trace=3, verbose=F)
}
names(FITS) <- names(DATA)

# save(FITS, file="data/HR_ctmm/FITS-SEASONS_all_skunks_error_2023-01-03.rda")
###########################################################
#calculate AKDEs with DATA and FITS files
load("data/HR_ctmm/FITS-SEASONS_all_skunks_error_2023-01-03.rda")

AKDES <- akde(DATA, FITS, weights=T) #must have verbose=F
meta(AKDES, col=c(COL,'black'), sort=TRUE) #km^2

#summarize info for hr and core areas
hr <- data.frame()
corearea <- data.frame()

for(i in c("SG-001","SG-001_2","SG-003","SG-005","SG-006","SG-008","SG-009","SG-010","SG-011","SG-013","SG-015"))
{
  print(i)
  animal <- grep(names(DATA), pattern=paste(i, "[.]", sep=""), value=T)
  length(animal)
  
  if(length(animal) == 2)
  {
    SEASON <- list(DATA[[animal[1]]], DATA[[animal[2]]])
    A <- list(AKDES[[animal[1]]], AKDES[[animal[2]]])
  }
  
  if(length(animal) == 3)
  {
    SEASON <- list(DATA[[animal[1]]], DATA[[animal[2]]], DATA[[animal[3]]])
    A <- list(AKDES[[animal[1]]], AKDES[[animal[2]]], AKDES[[animal[3]]])
  }
  
  # animal <- names(DATA)[[i]]
  # print(animal)
  COL <- rainbow(length(SEASON))
  
  # tiff(filename = paste("Figures/ctmm_seasonal/", i,"_hr.tiff", sep=""), height=8, width=8, units="in", res=330)
  plot(SEASON, UD=A, col=COL, lwd=2, col.grid=NA, col.level=COL, col.DF=NA, level.UD=0.95)
  # dev.off()
  
  # overlap(A)$CI[,,"est"]
  for(l in 1:length(animal))
  {
    hr <- rbind(hr, data.frame(animal[l], summary(A[[l]], level.UD=0.95)$CI, size="homerange"))
    corearea <- rbind(corearea, data.frame(animal[l], summary(A[[l]], level.UD=0.50)$CI, size="corearea"))
  }
}

hr <- separate(hr, col="animal.l.", into=c("animal","season"), sep="[.]", remove=F)
hr <- merge(hr, hr_ca_skunks[hr_ca_skunks$size=="homerange",c("animal","group")], by.x="animal", by.y="animal", all.x=T)
hr$season <- factor(hr$season, levels=c("spring","summer","fall"))
hr <- merge(hr, num.locs, by.x=c("animal","season"), by.y=c("animalid","season"), all.x=T)

corearea[grep(row.names(corearea), pattern="hectares"),]$low <- corearea[grep(row.names(corearea), pattern="hectares"),]$low/100
corearea[grep(row.names(corearea), pattern="hectares"),]$est <- corearea[grep(row.names(corearea), pattern="hectares"),]$est/100
corearea[grep(row.names(corearea), pattern="hectares"),]$high <- corearea[grep(row.names(corearea), pattern="hectares"),]$high/100
corearea <- separate(corearea, col="animal.l.", into=c("animal","season"), sep="[.]", remove=F)
corearea <- merge(corearea, hr_ca_skunks[hr_ca_skunks$size=="homerange",c("animal","group")], by.x="animal", by.y="animal", all.x=T)
corearea$season <- factor(corearea$season, levels=c("spring","summer","fall"))
corearea <- merge(corearea, num.locs, by.x=c("animal","season"), by.y=c("animalid","season"), all.x=T)

table(hr[hr$Freq > 6,c("season","group")])
#        group
# season F1 M1 M2
# spring  2  2  4
# summer  3  3  2
# fall    4  3  3

#compare changes in home range size by animal
#11 animals
hr.wide <- pivot_wider(hr[,c("animal", "season", "est", "group")], names_from = season, values_from=est)
hr.n.wide <- pivot_wider(hr[,c("animal", "season", "Freq", "group")], names_from = season, values_from=Freq)

hr <- merge(hr, skunk.info.unique, by.x="animal", by.y="Animal.ID", all.x=T)
hr$season <- factor(hr$season, levels=c("summer","fall","spring"))
hr <- hr[order(hr$season),]

p <- ggplot(hr, aes(x=season, y=est, col=animal, group=season)) + 
  geom_line(aes(group=animal), lty="dashed", size=1, position=position_dodge2(width=0.25)) + 
  geom_point(aes(size=mass), position=position_dodge2(width=0.25)) + 
  geom_errorbar(aes(ymin=low, ymax=high), linewidth=0.5, width=0.25, position=position_dodge2(width=1)) +
  theme_bw(base_size=20) + facet_wrap(~group, scales="free_x") + theme(panel.grid = element_blank()) + xlab("number of relocations") + ylab("home range size (km^2)") + coord_cartesian(ylim=c(0,75)) #by number of relocations
p
# ggsave(p, filename="Figures/ctmm_seasonal/seasonal_hr_changes_by_animal.tiff", dpi=400, height=7, width=12, units="in", compression="lzw")

########################################################
#HR overlap between seasons

o <- overlap(AKDES)
# write.csv(o$CI[,,"est"], file="data/HR_ctmm/overlap_seasonal_estimates_2023-01-03.csv")

o.est <- data.frame(o$CI[,,"est"])
o.est$animalid <- row.names(o.est)
o.est.long <- pivot_longer(o.est, !animalid, names_to="animalid2", values_to="overlap") #convert wide to long
o.est.long$animalid2 <- gsub(o.est.long$animalid2, pattern="SG.", replacement="SG-")

#same animal
o.same <- o.est.long[o.est.long$animalid != o.est.long$animalid2,]
o.same <- separate(o.same, col="animalid", remove=F, sep="[.]", into=c("a1","season1"))
o.same <- separate(o.same, col="animalid2", remove=F, sep="[.]", into=c("a2","season2"))

o.same <- o.same[o.same$a1 == o.same$a2,]
o.same$season <- paste(o.same$season1, "\n", o.same$season2, sep="")
o.same <- merge(o.same, hr_ca_skunks[hr_ca_skunks$size=="homerange",c("animal","group")], by.x="a1", by.y="animal", all.x=T)

p <- ggplot(o.same[o.same$season %in% c("spring\nsummer","spring\nfall","summer\nfall"),]) + 
  geom_path(aes(x=season, y=overlap*100, group=a1), position=position_dodge2(width=0.25)) + geom_point(aes(x=season, y=overlap*100, col=season), position=position_dodge2(width=0.25), size=3) + 
  ylab("percent overlap between seasons") + theme_bw(base_size=20) + facet_wrap(~group) + theme(legend.position = "none")
p

ddply(o.same, .(season), summarize, mean=mean(overlap), sd=sd(overlap), n=length(overlap), se=sd(overlap)/sqrt(length(overlap)))

3900/1824 #largest male average vs. largest male eastern spotted skunk home range: Lesmeister et al. 2009
#2.1
1630/1824
#0.89

1140/192 #largest female eastern spotted skunk home range: Lesmeister et al. 2009
#5.9

1630/61.1 #vs. MCP in dry season: Crooks and van vuren 1995
#26.7
1140/61.1
#18.66