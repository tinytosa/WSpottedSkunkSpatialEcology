#################################
#multi-season occupancy models to merge 2017 and 2018 data using package unmarked
#to determine best extinction model
#fixed year to reflect calendar year

#################################
#load packages
require(dplyr)
require(plyr)
require(unmarked)
require(AICcmodavg)

require(ggplot2)
require(ggpubr)

#################################
#load data
sitedata <- read.table("data/occupancy/carnivore_camera_allcovars.txt", sep=",") #includes siteinfo, buffer_logged, buffer_ogsi, landsat8, sentinel2, lidar, dem, gnn

sitedata$SiteName <- gsub(sitedata$Site, pattern="BP-", replacement="")
sitedata <- sitedata[!sitedata$Site %in% grep(sitedata$Site, pattern="SM-0[0-9]{1}[A|B|C]{1}", value = T),]

sitedata <- data.frame(sitedata, 'log_r' = log(sitedata$r.d + 1), 'log_w' = log(sitedata$w.d+1), 'log_s' = log(sitedata$s.d+1), 'e2' = sitedata$e^2) #log distance to road, distance to waterbody, distance to large stream
sitedata <- sitedata[order(sitedata$SiteName),]

sitedata[,grep(names(sitedata), pattern="percentage")] <- sitedata[,grep(names(sitedata), pattern="percentage")]/100

sitedata[is.na(sitedata$tpi.500m),]$tpi.500m <- mean(sitedata$tpi.500m, na.rm=T) #replace NA values with mean
sitedata[is.na(sitedata$tpi.1000m),]$tpi.1000m <- mean(sitedata$tpi.1000m, na.rm=T) #replace NA values with mean
sitedata[is.na(sitedata$Cover_4m_16m_all),]$Cover_4m_16m_all <- mean(sitedata$Cover_4m_16m_all, na.rm=T) #replace NA values with mean

#calculate mature (ogsi80 but not ogsi200)
sitedata$p_mature.r100 <- sitedata$percentage_inside.r100.o80 - sitedata$percentage_inside.r100.o200
sitedata$p_mature.r500 <- sitedata$percentage_inside.r500.o80 - sitedata$percentage_inside.r500.o200
sitedata$p_mature.r1000 <- sitedata$percentage_inside.r1000.o80 - sitedata$percentage_inside.r1000.o200
sitedata$p_mature.r5000 <- sitedata$percentage_inside.r5000.o80 - sitedata$percentage_inside.r5000.o200

##################
#detection histories

dhSpecies <- read.csv("data/occupancy/dh_WESTERN-SPOTTED-SKUNK_2022-10-24.csv")
row.names(dhSpecies) <- dhSpecies$SiteName

dhSpecies[!dhSpecies$SiteName %in% sitedata$SiteName,]$SiteName
dhSpecies <- dhSpecies[dhSpecies$SiteName %in% sitedata$SiteName,]

dhSpecies[dhSpecies > 0] <- 1 #convert to binary

##################
M <- nrow(dhSpecies) #112 sites
J <- 18 #number of secondary sample periods (weeks per season)
t <- 7 #number of primary sample periods (number of seasons)

#make 17 and 18 week seasons (really 18 week seasons)
#2017
dh2017.summer <- dhSpecies[,c(paste("X", 6:22, sep=""))] #6/3/2017 - 9/29/2017
dh2017.summer$XA <- NA

dh2017.fall <- dhSpecies[,c(paste("X", 23:39, sep=""))] #9/30/2017 - 1/26/2018
dh2017.fall$XB <- NA

dh2017.spring <- dhSpecies[,c(paste("X", 40:57, sep=""))] #1/27/2018 - 6/1/2018

#2018
dh2018.summer <- dhSpecies[,c(paste("X", 58:74, sep=""))] #6/2/2018 - 9/28/2018
dh2018.summer$XC <- NA

dh2018.fall <- dhSpecies[,c(paste("X", 75:91, sep=""))] #9/29/2018 - 1/25/2019
dh2018.fall$XD <- NA

dh2018.spring <- dhSpecies[,c(paste("X", 92:109, sep=""))] #1/26/2019 - 5/31/2019

#2019
dh2019.summer <- dhSpecies[,c(paste("X", 110:121, sep=""))] #6/1/2019 - 8/17/2019
dh2019.summer$XE <- NA
dh2019.summer$XF <- NA
dh2019.summer$XG <- NA
dh2019.summer$XH <- NA
dh2019.summer$XI <- NA
dh2019.summer$XJ <- NA

sp.matrix <- cbind(dh2017.summer, dh2017.fall, dh2017.spring, 
                   dh2018.summer, dh2018.fall, dh2018.spring, dh2019.summer)
dim(sp.matrix)

change <- data.frame(o2017.summer = rowSums(dh2017.summer, na.rm = T),
                     o2017.fall = rowSums(dh2017.fall, na.rm = T),
                     o2017.spring = rowSums(dh2017.spring, na.rm = T),
                     o2018.summer = rowSums(dh2018.summer, na.rm = T),
                     o2018.fall = rowSums(dh2018.fall, na.rm = T),
                     o2018.spring = rowSums(dh2018.spring, na.rm = T),
                     o2019.summer = rowSums(dh2019.summer, na.rm = T))
change[change > 0] <- 1

change$site <- row.names(change)
colSums(change[,-ncol(change)])

# o2017.summer   o2017.fall o2017.spring o2018.summer   o2018.fall o2018.spring o2019.summer 
#           31           39           28           44           57           41           20 
#        0.277        0.348        0.250        0.393        0.509        0.366        0.179 #if dividing by total sites for all seasons
#        0.574        0.709        0.509        0.393        0.509        0.366        0.179 #if dividing by number of sites that could have detected skunks

change.long <- tidyr::pivot_longer(change, cols=starts_with("o"), names_to="season", values_to="occupied")
change.long$season <- factor(change.long$season, levels=c("o2017.summer","o2017.fall","o2017.spring","o2018.summer","o2018.fall","o2018.spring","o2019.summer"))

p <- ggplot(change.long) + geom_tile(aes(x=season, y=site, fill=factor(occupied))) + scale_fill_manual(values=c("grey50","yellow")) + theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p
# ggsave(p, filename="Figures/multi-year_mult-season_occupancy_WSSK/occupancy_naive_byseason_2023-04-24.tiff", height=15, width=10, units="in", dpi=330, compression="lzw")

###########################
#observation covariates
#M*T*J rows
#112 * 18 * 7
#112 * 126

s.data <- read.table(file="data/occupancy/carnivore_camera_obs_covariates.txt", sep=",")
s.data <- s.data[s.data$week %in% gsub(colnames(sp.matrix), pattern="X", replacement=""),] #remove weeks not included in detection history

blankweeks <- data.frame(week=paste(LETTERS[1:10], sep=""), month=-1, num=-1, season="", weekofyear=-1, year=c(2016,2017,2017,2018,rep(2019,6)))

s.data$weekofyear <- 0
s.data[s.data$week < 36,]$weekofyear <- s.data[s.data$week < 36,]$week + 17 #calculate week of year
s.data[s.data$week >= 36,]$weekofyear <- s.data[s.data$week >= 36,]$week - 35 #calculate week of year
s.data[s.data$week >= 88,]$weekofyear <- s.data[s.data$week >= 88,]$week - 87 #calculate week of year
table(s.data$weekofyear)

s.data$year <- c(rep(2016, 30), rep(2017, 52), rep(2018, 34))
s.data <- rbind(s.data, blankweeks)
s.data$xweek <- paste("X", s.data$week, sep="")
s.data <- s.data[match(colnames(sp.matrix), s.data$xweek),] #reorder s.data to match order of sp.matrix

week.df <- data.frame(matrix(data=rep(as.character(s.data$weekofyear), times=nrow(sp.matrix)), ncol=length(s.data$week), byrow=T))
month.df <- data.frame(matrix(data=rep(as.character(s.data$month), times=nrow(sp.matrix)), ncol=length(s.data$month), byrow=T))
#season.df <- data.frame(matrix(data=rep(as.character(s.data$season), times=nrow(sp.matrix)), ncol=length(s.data$season), byrow=T)) #need equal sized number of weeks for each season
calyear.df <- data.frame(matrix(data=rep(as.character(s.data$year), times=nrow(sp.matrix)), ncol=length(s.data$year), byrow=T))


#time (weeks) since bait
bait.blank <- data.frame(matrix(rep(-1, 10*M), ncol=10))
names(bait.blank) <- paste("X", LETTERS[1:10], sep="")

bait.df <- read.csv("data/occupancy/carnivore_camera_weekssincebait.csv")
bait.df <- bait.df[bait.df$site %in% sitedata$SiteName,] #remove cameras not included in detection history
bait.df <- bait.df[,names(bait.df) %in% colnames(sp.matrix)] #remove weeks not included in detection history

bait.df <- cbind(bait.df, bait.blank) #add in bait.blank
bait.df <- bait.df[,colnames(sp.matrix)] #reorder bait.df
dim(bait.df)

###########################
#yearlySiteCovs
#should be a data frame with M*t rows
year <- matrix(c('2017','2017','2017','2018','2018','2018','2019'), M, t, byrow=T)
season <- matrix(c('summer','fall','spring','summer','fall','spring','summer'), M, t, byrow=T)

##########################
#create unmarked frame
umf <- unmarkedMultFrame(y=sp.matrix, siteCovs=sitedata, 
                         yearlySiteCovs=list(year=year, season=season),
                         obsCovs=list(week=week.df, month=month.df, calyear=calyear.df, bait=bait.df), numPrimary=t) #has unscaled siteCovs
head(umf)
plot(umf)
summary(umf)

##########################
#run occupancy models
# ~psi, ~gamma, ~ epsilon, ~p, data
#psi = occupancy during first season
#gamma = colonization rate
#epsilon = extinction rate
#p = detection probability
##########################

p.global <- colext(~1, ~1, ~1, ~season + year + bait, data=umf) #can't use week or month because estimates are weird
psi.greduced <- colext(~scale(p_mature.r5000) + scale(B4_20180818) + scale(tpi.1000m), ~1, ~1, ~1, data=umf)
gamma.greduced <- colext(~1, ~season + scale(r.d) + scale(B6_20180818) + scale(tri) + scale(abam_ba_2017), ~1,  ~1, data=umf)
e.greduced <- colext(~1, ~1, ~season + year + scale(p_mature.r5000) + scale(YrsSinceDist) + log(YrsSinceDist + 0.001) + 
                       scale(tpi.500m) + scale(e) + scale(e^2) + scale(abam_ba_2017), ~1, data=umf)

global <- colext(~scale(p_mature.r5000) + scale(B4_20180818) + scale(tpi.1000m),
                 # ~scale(r.d) +
                 ~season + scale(B6_20180818) + scale(tri) + scale(abam_ba_2017), 
                 ~season + year + scale(p_mature.r5000) + scale(YrsSinceDist) + log(YrsSinceDist + 0.001) + scale(tpi.500m) + scale(e) + scale(e^2) + scale(abam_ba_2017),
                 ~season + year + bait, data=umf)

summary(global)

mod <- global
coefs <- coef(mod) %>% tibble::enframe()
se <- SE(mod) %>% tibble::enframe()

ci.psi <- confint(mod, type="psi") %>% tibble::enframe()
ci.col <- confint(mod, type="col") %>% tibble::enframe()
ci.ext <- confint(mod, type="ext") %>% tibble::enframe()
ci.det <- confint(mod, type="det") %>% tibble::enframe()

ci <- rbind(ci.psi, ci.col, ci.ext, ci.det)
ci <- data.frame(name=ci$name, lower=ci$value[,1], upper=ci$value[,2])

coefs <- merge(coefs, se, by="name")
coefs <- merge(coefs, ci, by="name")
names(coefs) <- c("par", "est", "se", "lower","upper")

coefs$p <- gsub(substr(coefs$par, 0,3), pattern="p\\([A-z]", replacement="det")

coefs$name <- gsub(coefs$par, pattern="ext\\(", replacement="")
coefs$name <- gsub(coefs$name, pattern="scale\\(", replacement="")
coefs$name <- gsub(coefs$name, pattern="\\)", replacement="")
coefs$name <- gsub(coefs$name, pattern="_2017", replacement="")
coefs$name <- gsub(coefs$name, pattern=" \\+ 0.001", replacement="\\)")
coefs$name <- gsub(coefs$name, pattern="\\_GE_3", replacement="", ignore.case = T)
coefs$name <- gsub(coefs$name, pattern="_20180818", replacement="")
coefs$name <- gsub(coefs$name, pattern="percentage_inside[.]", replacement="p.logged.")
coefs$name <- gsub(coefs$name, pattern="col\\(", replacement="")
coefs$name <- gsub(coefs$name, pattern="p\\(", replacement="")
coefs$name <- gsub(coefs$name, pattern="psi\\(", replacement="")
coefs$name <- gsub(coefs$name, pattern="_ba", replacement="")
coefs$name <- gsub(coefs$name, pattern="r.d", replacement="dist.road")
coefs$name <- gsub(coefs$name, pattern="r5000", replacement="5km")
coefs$name <- gsub(coefs$name, pattern="r1000", replacement="1km")
coefs$name <- gsub(coefs$name, pattern="1000m", replacement="1km")
coefs$name <- gsub(coefs$name, pattern="season", replacement="")
coefs$name <- gsub(coefs$name, pattern="year", replacement="")
coefs$name <- gsub(coefs$name, pattern="^e$", replacement="elevation")
coefs$name <- gsub(coefs$name, pattern="tpi.", replacement="Topo pos.")
coefs <- coefs[coefs$name != "Int",]

coefs$name <- factor(toupper(coefs$name), levels=toupper(c("Int","ogsi","log(ogsi)","YrsSinceDist","log(YrsSinceDist)","Dist.road","p.logged.5km","p.mature.100m","p_mature.5km","p_oldgrowth.0.5km",
                                                           "B2","B3","B4","B5","B6","B7",
                                                           "Dist.stream","Topo pos.0.5km","Topo pos.500m","Topo pos.1km",
                                                           "Cover","Canopy cover","tphc","Conifer density","slope","cancov_con","tri","Rough",
                                                           "elevation","e2","elevation2","e^2",
                                                           "ABAM","abam_ba","Abies amabilis","ACMA","PSME","Pseudotsuga menzeseii","bah",
                                                           "2018","2019","spring","summer","bait")))
coefs$category <- c("Thermal","Resource","Predation", rep("Temporal", 2),
                    "Disturbance","Thermal","Thermal","Thermal","Disturbance","Resource","Disturbance",
                    rep("Temporal", 8),
                    "Resource","Disturbance","Resource")
coefs$p <- factor(coefs$p, levels=c("det","psi","col","ext"))

# write.csv(coefs, file="data/occupancy/occupancy_coefficients.csv", row.names=F)

c.plot <- ggplot(data=coefs, aes(x=name, y=est, group=category, shape=category)) +
  geom_point(size=3, position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin=est-se, ymax=est+se), width=0, lwd=1, position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept=0), lty="dotted") +
  xlab("parameter") + ylab("estimate") +
  coord_flip() +
  theme_bw(base_size=16) + theme(panel.grid=element_blank()) + facet_wrap(.~p, scales="free", ncol=1)
c.plot

# ggsave(c.plot, filename=paste("Figures/multi-year_mult-season_occupancy_WSSK/final_betacoefficients_all.tiff", sep=""), height=10, width=8, units="in", dpi=300, compression="lzw")

###########################
#plot detection probabilities for final model
# season + year + bait

nd.og <- data.frame(bait=0, year='2017', season='summer')
nd2 <- data.frame(nd.og[,names(nd.og) != "bait"], bait=0:max(bait.df))
marginal2 <- predict(global, type="det", newdata=nd2)
d2 <- ggplot(marginal2, aes(x=nd2$bait, y=Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
  ylab("Detection probability") +
  geom_line(lwd=1) +theme_bw(base_size=30) + xlab("weeks since bait") + theme(panel.grid = element_blank())
d2
# ggsave(d2, filename="Figures/multi-year_mult-season_occupancy_WSSK/det_global_bait.tiff", height=8, width=12, units="in", dpi=300, compression="lzw")

#put categorical variables in the same figure
nd <- data.frame(bait=0, year=c('2017','2017','2017','2018','2018','2018','2019'), season=c('summer','fall','spring','summer','fall','spring','summer'))
marginal <- predict(global, type="det", newdata=nd)
marginal <- cbind(marginal, nd)
marginal$calyear <- c('2017','2017','2018','2018','2018','2019','2019')
marginal$season <- factor(marginal$season, levels=c("summer","fall","spring"))
d <- ggplot(marginal, aes(x=season, y=Predicted, col=year, shape=year)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), size=2, width=0.1, position = position_dodge(width=0.5)) +
  ylim(c(0,NA)) +
  xlab("season") + ylab("Detection probability") +
  scale_shape_manual(values=c(16,17,15), name="skunk year") +
  scale_color_manual(values=c("#334E58","#80A4ED","#84596B"), name="skunk year") +
  geom_point(size=10, position = position_dodge(width=0.5)) + theme_bw(base_size=30) + theme(panel.grid = element_blank(), legend.position = c(0.85, 0.2))
d
# ggsave(d, filename="Figures/multi-year_mult-season_occupancy_WSSK/det_global_temporal.tiff", height=8, width=12, units="in", dpi=300, compression="lzw")

marginal %>% mutate_if(is.numeric, round, 3)

#   Predicted    SE lower upper bait year season calyear
# 1     0.165 0.014 0.140 0.195    0 2017 summer    2017
# 2     0.356 0.018 0.322 0.391    0 2017   fall    2017

# 3     0.241 0.019 0.205 0.280    0 2017 spring    2018
# 4     0.203 0.015 0.175 0.234    0 2018 summer    2018
# 5     0.415 0.017 0.383 0.448    0 2018   fall    2018

# 6     0.289 0.021 0.250 0.332    0 2018 spring    2019
# 7     0.222 0.036 0.159 0.301    0 2019 summer    2019

# ggsave(ggarrange(d2, d, nrow=1, align="h"), filename=paste("Figures/multi-year_mult-season_occupancy_WSSK/final_detection_monthbait.tiff", sep=""), height=8, width=16, unit="in", dpi=300, compression="lzw")


###############################
#marginal plots for each variable

#plot non-scaled variables
#create new data frame with means of all covariates

########
#plot marginals for 
nd.mean <- sitedata %>% mutate_if(is.numeric, mean)
nd.mean <- nd.mean[1,] #only take first row since all rows are the same
nd.mean$season <- "fall"
nd.mean$year <- "2017"

ylabel <- data.frame(type=c("det","psi","col","ext"), label=c("detection","occupancy","colonization","extinction"))
#
plot.marginal <- function(par, par.name, type)
{
  nd <- data.frame(nd.mean[,names(nd.mean) != par], par=seq(min(sitedata[,par]), max(sitedata[,par]), length.out=1000))
  names(nd) <- gsub(names(nd), pattern="par", replacement=par)
  marginal <- predict(global, type=type, newdata=nd)
  
  o <- ggplot(marginal, aes(x=nd[,par], y=Predicted)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
    geom_rug(data=sitedata, aes(x=sitedata[,par], y=0)) +
    ylim(0,1) +
    xlab(par.name) + ylab(paste("Predicted probability of", ylabel[ylabel$type == type,]$label)) + 
    geom_line() + theme_bw(base_size=20) + theme(panel.grid=element_blank())
  print(o)
  # ggsave(o, filename=paste("Figures/multi-year_mult-season_occupancy_WSSK/final_",type,"_", par,".tiff", sep=""), height=6, width=6, unit="in", dpi=300, compression="lzw")
}


plot.marginal("B4_20180818", "LANDSAT B4", type="psi") # (-)
plot.marginal("tpi.1000m", "Topographic position index within 1km", type="psi") # (-)
plot.marginal("p_mature.r5000", "percent mature within 5km", type="psi") # (+)

plot.marginal("B6_20180818", "Landsat B6", type="col")
plot.marginal("tri", "Roughness", type="col")
plot.marginal("abam_ba_2017", "Abies amabilis basal area", type="col")

plot.marginal("abam_ba_2017","Abies amabilis basal area", type="ext")
plot.marginal("tpi.500m","Topographic position index (0.5km)", type="ext")
plot.marginal("YrsSinceDist","Years since disturbance", type="ext")
plot.marginal("e", "Elevation (m)", type="ext")
plot.marginal("p_mature.r5000", "Percent mature within 5km", type="ext")

#new function in unmarked! only plot in base R, not ggplot
plotEffects(mod, "psi", "B4_20180818")

########
#plot marginals for factor variables
nd.mean <- sitedata %>% mutate_if(is.numeric, mean)
nd.mean <- nd.mean[1,] #only take first row since all rows are the same

nd <- nd.mean
nd <- data.frame(nd, season=c("summer","fall","spring"))
nd$season <- factor(nd$season, levels=c("summer","fall","spring"))
marginal <- predict(global, type="col", newdata=nd)

o <- ggplot(marginal, aes(x=nd[,"season"], y=Predicted)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), size=1, width=0.1) +
  ylim(0,1) +
  xlab("season") + ylab("Predicted probability of colonization") + scale_x_discrete(breaks=c("summer","fall","spring"), labels=c("summer-fall","fall-spring","spring-summer")) +
  geom_point(size=5) + theme_bw(base_size=20) + theme(panel.grid=element_blank())
print(o)
# ggsave(o, filename=paste("Figures/multi-year_mult-season_occupancy_WSSK/final_","col","_", "season",".tiff", sep=""), height=6, width=6, unit="in", dpi=300, compression="lzw")

#year
nd <- nd.mean
nd <- data.frame(nd, year=c('2017','2017','2017','2018','2018','2018'), season=c('summer','fall','spring','summer','fall','spring'))
marginal <- predict(global, type="ext", newdata=nd)
marginal <- data.frame(marginal, year=nd$year, season=nd$season)
marginal$season <- factor(marginal$season, levels=c("summer","fall","spring"))

o <- ggplot(marginal, aes(x=season, y=Predicted, col=year, shape=year)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), size=1, width=0.1, position = position_dodge(width=0.5)) +
  scale_color_manual(values=c("#334E58","#80A4ED","#84596B"), name="skunk year") +
  scale_shape_manual(values=c(16,17,15), name="skunk year") +
  ylim(0,1) +
  xlab("season") + ylab("Predicted probability of extinction") + scale_x_discrete(breaks=c("summer","fall","spring"), labels=c("summer-fall","fall-spring","spring-summer")) +
  geom_point(size=7.5, position = position_dodge(width=0.5)) + theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position = c(0.8, 0.8))
print(o)
# ggsave(o, filename=paste("Figures/multi-year_mult-season_occupancy_WSSK/final_","ext","_", "season-year",".tiff", sep=""), height=6, width=6, unit="in", dpi=300, compression="lzw")


###########################
################
# Mackenzie-Bailey GOF test
# from https://jamesepaterson.github.io/jamespatersonblog/2021-01-01_dynamicoccupancy.html
# Simulate capture history data (if model correct). Compare obs X2 to sim X2
# Must simulate 1000-5000 times for good distribution estimates
# Likely to take a Very Long Time on large data sets
# with nsim = 10, takes my Mackbook ~ 60 seconds
mb.boot <- AICcmodavg::mb.gof.test(global, nsim = 1000) # Must be much higher than five to be useful
mb.boot

# Goodness-of-fit for dynamic occupancy model
# 
# Number of seasons:  7 
# 
# Chi-square statistic:
#   Season 1     Season 2     Season 3     Season 4     Season 5     Season 6     Season 7 
# 43364329.422  6650121.807  1859665.959 24728076.826  1181136.112 17132951.361     3953.516 
# 
# Total chi-square = 94920235 
# Number of bootstrap samples = 1000
# P-value = 0
# 
# Quantiles of bootstrapped statistics:
#   0%     25%     50%     75%    100% 
# 3.0e+05 5.5e+05 6.9e+05 1.0e+06 8.6e+07 
# 
# Estimate of c-hat = 79.2 

# from unmarked vignette: https://rdrr.io/cran/unmarked/f/vignettes/unmarked.Rmd
# The parametric bootstrap can be used to check the adequacy of model fit. Here we use a $\chi^2$ statistic appropriate for binary data.
chisq <- function(fm) {
  umf <- fm@data
  y <- umf@y
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

(pb <- parboot(global, statistic=chisq, nsim=1000, parallel=FALSE))

#predicted occupancy per site per season
global@smoothed[1,,] #1 for unoccupied
global@smoothed[2,,] #2 for occupied

#mean predicted occupancy among sites
global@smoothed.mean

#                    1         2         3         4         5         6         7
# unoccupied 0.4983108 0.3626121 0.4521918 0.6067561 0.4908031 0.6053123 0.8008715
# occupied   0.5016892 0.6373879 0.5478082 0.3932439 0.5091969 0.3946877 0.1991285

global@projected.mean
#                    1         2         3         4         5         6         7
# unoccupied 0.5003834 0.3915741 0.4614082 0.6171883 0.4782339 0.6026794 0.7849925
# occupied   0.4996166 0.6084259 0.5385918 0.3828117 0.5217661 0.3973206 0.2150075

#non-parametric bootstrap to obtain standard errors of smoothed estimates of occupancy probability during each season
global <- nonparboot(global, B=1000)
smooth <- cbind(smoothed=smoothed(global)[2,], SE=global@smoothed.mean.bsse[2,])
smooth <- data.frame(smooth)
smooth$season <- factor(c('summer\n2017','fall\n2017','spring\n2018','summer\n2018','fall\n2018','spring\n2019','summer\n2019'),
                        levels=c('summer\n2017','fall\n2017','spring\n2018','summer\n2018','fall\n2018','spring\n2019','summer\n2019'))
smooth$num <- as.numeric(row.names(smooth))
smooth

#    smoothed         SE       season num
# 1 0.5016892 0.08785094 summer\n2017   1
# 2 0.6373879 0.05312252   fall\n2017   2
# 3 0.5478082 0.05574598 spring\n2018   3
# 4 0.3932439 0.04850887 summer\n2018   4
# 5 0.5091969 0.04775449   fall\n2018   5
# 6 0.3946877 0.04638620 spring\n2019   6
# 7 0.1991285 0.04265960 summer\n2019   7

o.plot <- ggplot(smooth, aes(x=num, y=smoothed)) +
  geom_smooth(aes(x=num, y=smoothed), method="lm", col="black", lty="dashed") +
  geom_smooth(aes(x=num, y=smoothed), method="lm", col=NA, lty="dashed", fill=NA) +
  geom_point(size=3) + 
  geom_errorbar(aes(ymin=smoothed-SE, ymax=smoothed+SE), width=0.25, linewidth=1) + 
  scale_x_continuous(breaks=1:7, labels=smooth$season) +
  xlab("season") + ylab("predicted occupancy") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank())
o.plot

# ggsave(o.plot, filename="Figures/multi-year_mult-season_occupancy_WSSK/predicted_occupancy_across_sites_smoothed.tiff", #with trend line
#        # ggsave(o.plot, filename="Figures/multi-year_mult-season_occupancy_WSSK/predicted_occupancy_across_sites.tiff", #without trendline
#        height=6, width=8, units="in", dpi=300, compression="lzw")

########################
#fit linear regression to smooth
smooth <- read.table(text="
smoothed SE season num
0.5016892 0.08785094 summer_2017 1
0.6373879 0.05312252 fall_2017 2
0.5478082 0.05574598 spring_2018 3
0.3932439 0.04850887 summer_2018 4
0.5091969 0.04775449 fall_2018 5
0.3946877 0.04638620 spring_2019 6
0.1991285 0.04265960 summer_2019 7", sep=" ", header=T)

smooth$num <- as.numeric(row.names(smooth))
summary(lm(smoothed ~ num, data=smooth))

smooth$season <- gsub(smooth$season, pattern="_", replacement="\n")

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)  0.65926    0.08169   8.070 0.000473 ***
#   num         -0.05113    0.01827  -2.799 0.038027 *  

o.plot <- ggplot(smooth, aes(x=num, y=smoothed)) +
  geom_smooth(aes(x=num, y=smoothed), method="lm", col="black", lty="dashed") +
  geom_smooth(aes(x=num, y=smoothed), method="lm", col=NA, lty="dashed", fill=NA) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=smoothed-SE, ymax=smoothed+SE), width=0.25, linewidth=1) +
  scale_x_continuous(breaks=1:7, labels=smooth$season) +
  xlab("season") + ylab("predicted occupancy") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank())
o.plot

# ggsave(o.plot, filename="Figures/mult-season-season_occupancy_WSSK/predicted_occupancy_across_sites.tiff", #without trendline
#        height=6, width=8, units="in", dpi=300, compression="lzw")

