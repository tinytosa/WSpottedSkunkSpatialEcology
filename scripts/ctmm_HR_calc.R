#######################
#skunk data with ctmm
#######################

#load packages
require(plyr)
require(tidyr)
require(dplyr)
require(ctmm)
require(rgdal)

require(ggplot2)
require(ggpubr)

###########
#load data


##########
#VHF data
#gave rest sites COV.x.x <- 1, COV.y.y <- 1, COV.x.y <- 0

nad83z10 <- "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

sk.summary <- read.csv(file="data/sk_summary.csv", header=T)
sk.locs <- read.table(file="data/HR_ctmm/as_telemetry_obj_vhf_2022-11-17.txt", sep=",", header=T)
table(sk.locs$individual.local.identifier)

DATA <- as.telemetry(sk.locs, projection=CRS(nad83z10))

uere(DATA) <- 1 # you need this in newer versions of the package!

#read in error ellipse information for each skunk
for(i in 1:length(DATA))
{
  animal <- names(DATA)[[i]]
  print(animal)
  
  DATA[[i]]$COV.x.x <- sk.locs[sk.locs$individual.local.identifier == animal,]$COV.x.x
  DATA[[i]]$COV.y.y <- sk.locs[sk.locs$individual.local.identifier == animal,]$COV.y.y
  DATA[[i]]$COV.x.y <- sk.locs[sk.locs$individual.local.identifier == animal,]$COV.x.y
  DATA[[i]]$VAR.xy <- (DATA[[i]]$COV.x.x + DATA[[i]]$COV.y.y)/2
}

#plot vhf locations with error ellipses
COL <- rainbow(length(DATA))
plot(DATA, col=COL, lwd=2)

###########
#GPS data

gps.locs1 <- read.table("data/HR_ctmm/GPS_SG-011.txt", sep=",", header=T)
gps.locs2 <- read.table("data/HR_ctmm/GPS_SG-019.txt", sep=",", header=T)

GPS1 <- as.telemetry(gps.locs1, projection=CRS(nad83z10))
GPS2 <- as.telemetry(gps.locs2, projection=CRS(nad83z10))

#calibration of data with a prior
summary(uere(GPS1))
uere(GPS1) <- 20
UERE <- uere(GPS1)
summary(UERE)
UERE$DOF[] <- 2
summary(UERE)
uere(GPS1) <- UERE

#calibration of data with a prior
summary(uere(GPS2))
uere(GPS2) <- 10
UERE <- uere(GPS2)
summary(UERE)
UERE$DOF[] <- 2
summary(UERE)
uere(GPS2) <- UERE

###############
#combine vhf with gps locations
DATA[['SG-011']] <- tbind(DATA[['SG-011']], GPS1)
DATA[['SG-019']] <- tbind(DATA[['SG-019']], GPS2)

plot(DATA, col=COL, lwd=2)
table(sk.locs$individual.local.identifier)

# save(DATA, file="data/HR_ctmm/DATA_all_skunks_error_2022-11-22.rda")

FITS <- list()
for(i in 1:length(DATA))
  # for(i in 1:3)
{
  animal <- names(DATA)[[i]]
  print(animal)
  
  # tiff(filename = paste("Figures/ctmm_points/plot_", animal,".tiff", sep=""), height=8, width=8, units="in", res=300)
  COL <- color(DATA[[i]], by='time')
  plot(DATA[[i]], col=COL, lwd=2)
  # dev.off()
  
  GUESS <- ctmm.guess(DATA[[i]], interactive=F, CTMM = ctmm(error=T)) #CTMM-ctmm(error=T) include this 
  FITS[[i]] <- ctmm.select(DATA[[i]], GUESS, trace=3, verbose=F)
  #add plot and save here to look at post hoc
  # SVF <- variogram(DATA[[i]])
  # tiff(filename = paste("Figures/ctmm_checks/SVF_", animal,".tiff", sep=""), height=8, width=8, units="in", res=300)
  #   plot(SVF, FITS[[i]][[1]])
  # dev.off()
  
  # RES2 <- residuals(DATA[[i]], FITS[[i]][[1]])
  # ACF2 <- correlogram(RES2, res=10)
  # 
  # tiff(filename = paste("Figures/ctmm_checks/ACF_", animal,".tiff", sep=""), height=8, width=8, units="in", res=300)
  #   plot(ACF2)
  # dev.off()
}
names(FITS) <- names(DATA)

# save(FITS, file="data/HR_ctmm/FITS_all_skunks_error_2022-11-16.rda")

#############
load("data/HR_ctmm/DATA_all_skunks_error_2022-11-22.rda")
COL <- rainbow(length(DATA))
load("data/HR_ctmm/FITS_all_skunks_error_2022-11-16.rda")

AKDES <- akde(DATA, FITS, weights=T) #must have verbose=F

# plot AKDEs
# tiff(filename="Figures/homeranges_ctmm_with_error_2022-11-16.tiff", height=8, width=8, units="in", res=300, compression="lzw")
COL <- color(AKDES, by='individual')
plot(AKDES, col.DF=COL, col.level=COL, col.grid=NA, level=NA)
# dev.off()

############
#95% isopleth home range information
# meta-analysis of skunk home-range areas
# tiff(filename="Figures/homeranges_ctmm_est_with_error_2022-11-16.tiff", height=8, width=8, units="in", res=300, compression="lzw")
meta(AKDES, col=c(COL,'black'), sort=TRUE) #km^2
# dev.off()

#95% UD areas
# ΔAICc
# inverse-Gaussian   0.0000
# Dirac-δ           361.6763
# 
#                    low        est      high
# mean (km²)   14.258607 19.9165839 26.975123
# CoV² (RVAR)   0.223124  0.5457124  1.010060
# CoV  (RSTD)   0.480473  0.7514112  1.022279

#####
#plot each home range akde with points, save core area and home range area information
#given in square kilometers
hr <- data.frame()
corearea <- data.frame()
for(i in 1:length(DATA))
{
  animal <- names(DATA)[[i]]
  print(animal)
  COL <- color(DATA[[i]], by='time')
  
  # tiff(filename = paste("Figures/ctmm_hr/", animal,"_hr.tiff", sep=""), height=8, width=8, units="in", res=300)
  #   plot(DATA[[i]], UD=AKDES[[i]], col=COL, lwd=2, col.grid=NA, level.UD=0.95)
  # dev.off()
  # 
  # tiff(filename = paste("Figures/ctmm_corearea/", animal,"_corearea.tiff", sep=""), height=8, width=8, units="in", res=300)
  #   plot(DATA[[i]], UD=AKDES[[i]], col=COL, lwd=2, col.grid=NA, level.UD=0.5)
  # dev.off()
  
  hr <- rbind(hr, data.frame(animal, summary(AKDES[[i]], level.UD=0.95)$CI, size="homerange"))
  corearea <- rbind(corearea, data.frame(animal, summary(AKDES[[i]], level.UD=0.50)$CI, size="corearea"))
}

#############
#cluster males based on separation
set.seed(1)
cluster <- kmeans(hr$est, centers=2, iter.max=100, nstart=2) #2 groups #this looks best
# cluster <- kmeans(hr$est, centers=3, iter.max=100, nstart=2) #3 groups
# cluster <- kmeans(hr$est, centers=4, iter.max=100, nstart=2) #4 groups
# cluster <- kmeans(hr$est, centers=5, iter.max=100, nstart=2) #5 groups

cluster.info <- data.frame(cluster=cluster$cluster, animal=hr$animal)

m1 <- hr[cluster$cluster == 1,]$animal #smaller male home ranges
# "SG-002" "SG-006" "SG-011" "SG-013" "SG-016" "SG-021" "SG-023" "SG-025"
m2 <- hr[cluster$cluster == 2,]$animal #larger male home ranges
# "SG-001"   "SG-001_2" "SG-003"   "SG-004"   "SG-005"   "SG-007"   "SG-008"   "SG-009"   "SG-010"   "SG-012"   "SG-014"   "SG-015"   "SG-017"   "SG-018"   "SG-019"   "SG-020"   "SG-022" 
f <- c("SG-005", "SG-007", "SG-008", "SG-009", "SG-015", "SG-016", "SG-017","SG-018","SG-020") #females

a <- rbind(hr, corearea)

a$sex <- "M"
a[a$animal %in% f,]$sex <- "F"

a[grep(row.names(a), pattern="hectares"),]$low <- a[grep(row.names(a), pattern="hectares"),]$low/100
a[grep(row.names(a), pattern="hectares"),]$est <- a[grep(row.names(a), pattern="hectares"),]$est/100
a[grep(row.names(a), pattern="hectares"),]$high <- a[grep(row.names(a), pattern="hectares"),]$high/100

a <- merge(a, sk.summary, by.x="animal", by.y="AnimalID", all.x=T)
a <- merge(a, cluster.info, by="animal", all.x=T)
a$cluster <- factor(a$cluster)
a$group <- paste(a$sex, a$cluster, sep="")
a[a$group == "F2",]$group <- "F1"
a$title <- paste(a$size, a$group)

# write.table(a, file="data/HR_ctmm/hr_ca_estimates.txt", sep=",", row.names=F)

########################################################################

######################
a <- read.table("data/HR_ctmm/hr_ca_estimates.txt", sep=",", header=T)
a$cluster <- factor(a$cluster)

#compare home range and core area size across groups
a[a$size == "homerange",] %>% group_by(group) %>% dplyr::summarize(mean = mean(est, rm.na=T), n=length(est), sd=sd(est), se = sd(est)/sqrt(length(est)))
#   group  mean     n    sd    se
#   <chr> <dbl> <int> <dbl> <dbl>
# 1 F1     11.4     9  8.14  2.71
# 2 M1     39.0     7 16.6   6.27
# 3 M2     16.3     9  5.15  1.72
summary(aov(data=a[a$size == "homerange",], est ~ group))
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# group        2   3301  1650.7   15.17 7.23e-05 ***
# Residuals   22   2394   108.8                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov(data=a[a$size == "homerange",], est ~ group))
# $group
#            diff        lwr       upr     p adj
# M1-F1  27.61817  14.412921 40.823409 0.0000819
# M2-F1   4.88964  -7.462734 17.242015 0.5880346
# M2-M1 -22.72852 -35.933769 -9.523281 0.0007695

###
summary(aov(data=a[a$size == "corearea",], est ~ group))
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# group        2  190.5   95.27   11.98 0.000302 ***
# Residuals   22  174.9    7.95                   
TukeyHSD(aov(data=a[a$size == "corearea",], est ~ group))
#             diff       lwr       upr     p adj
# M1-F1  6.5566471  2.987306 10.125988 0.0003811
# M2-F1  0.9468259 -2.391987  4.285639 0.7588287
# M2-M1 -5.6098211 -9.179162 -2.040480 0.0019016


#########################################
#figures
#########################################
base <- ggplot(a, aes(x=est, y=animal)) +
  geom_point(size=1) + geom_errorbar(aes(xmin=low, xmax=high), lwd=1, width=0.5) +
  xlab("area (km^2)") + ylab("") + 
  theme_bw(base_size = 20)

base + geom_point(size=3, aes(col=cluster, shape=sex)) + facet_wrap(~size + sex, scale="free", nrow=2)

p.ests <- base + geom_point(size=3) +
  geom_errorbar(aes(xmin=low, xmax=high), lwd=1, width=0.5) +
  facet_wrap(~title, scale="free", nrow=3, dir = "v")
p.ests
# ggsave(p.ests, filename="Figures/ctmm_ests_2022-11-18.tiff", height=15, width=15, units="in", compression="lzw", dpi=330)

ggplot(a, aes(x=est, y=animal, col=cluster, shape=sex), position=position_dodge2(width = 2)) + 
  geom_point(size=5) + 
  geom_errorbar(aes(xmin=low, xmax=high), lwd=1, width=0.5) + 
  theme_bw(base_size = 20) + facet_wrap(~cluster + size, scales="free", nrow=3, ncol=2)

p.locs <- ggplot(a, aes(x=total, y=est, col=group, shape=sex)) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=low, ymax=high), lwd=0.5, width=10) +
  xlab("number of relocations") + ylab("area estimate (km^2)") +
  theme_bw(base_size = 20) + facet_wrap(~size, scale="free_y") + scale_y_log10()
p.locs
# ggsave(p.locs, filename = "Figures/ctmm_locs_vs_area.tiff", height=4, width=10, units="in", compression="lzw", dpi=330)
