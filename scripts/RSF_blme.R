##########
#skunk rsf at home range level
#univariate models


#try with blme since glmer yields var 0 & isSingular=T

#https://ms.mcmaster.ca/~bolker/R/misc/foxchapter/bolker_chap.html
# use bglmer from the blme package, which sets a weak prior on the variance to push it away from zero:

###########
#load packages
require(blme) #for bglmer
require(bbmle) #for AICtab

require(AICcmodavg)
require(ggplot2)
require(ggpubr)
require(lattice)
require(ggh4x)
require(ggeffects)

###########
#load data
#skunk data
skunk.info <- read.table("data/HR_ctmm/hr_ca_estimates.txt", sep=",", header=T)

#covariates
covars <- read.table("Data/RSF/covar_rsf_pts.txt", sep=",")
covars$elev <- covars$elev*0.3048 #convert from ft to m
covars.lidar <- read.table("Data/RSF/covar_rsf_lidar.txt", sep=",")
covars.landsat <- read.table("Data/RSF/covars_rsf_landsat8.txt", sep=",")
covars.sentinel <- read.table("Data/RSF/covars_rsf_sentinel2.txt", sep=",")
covars.distance <- read.table("Data/RSF/covars_rsf_distanceto.txt", sep=",")
covars.buffer <- read.table("Data/RSF/covar_rsf_buffer_logged.txt", sep=",")
covars.buffer.ogsi <- read.table("Data/RSF/covar_rsf_buffer_ogsi.txt", sep=",")

#calculate amount o80 but not o200 for each buffer size
covars.buffer.ogsi$p_mature.r100 <- covars.buffer.ogsi$percentage_inside.r100.o80 - covars.buffer.ogsi$percentage_inside.r100.o200
covars.buffer.ogsi$p_mature.r500 <- covars.buffer.ogsi$percentage_inside.r500.o80 - covars.buffer.ogsi$percentage_inside.r500.o200
covars.buffer.ogsi$p_mature.r1000 <- covars.buffer.ogsi$percentage_inside.r1000.o80 - covars.buffer.ogsi$percentage_inside.r1000.o200
covars.buffer.ogsi$p_mature.r5000 <- covars.buffer.ogsi$percentage_inside.r5000.o80 - covars.buffer.ogsi$percentage_inside.r5000.o200

remove <- c("AnimlID","UTM_E","UTM_N","type","used","pt")

#combine all covars
data <- cbind(covars,
              covars.lidar[,!names(covars.lidar) %in% remove],
              covars.landsat[,!names(covars.landsat) %in% remove],
              covars.sentinel[,!names(covars.sentinel) %in% remove],
              covars.distance[,!names(covars.distance) %in% remove],
              covars.buffer[,!names(covars.buffer) %in% remove],
              covars.buffer.ogsi[,!names(covars.buffer.ogsi) %in% remove]
)
dim(data)
#39832 points x 174 covariates
#39832 points x 178 covariates (with mature, not old)

#fill in missing data with mean values
data[is.na(data$tpi.1000m),]$tpi.1000m <- mean(data$tpi.1000m, na.rm=T)
data[is.na(data$tpi.500m),]$tpi.500m <- mean(data$tpi.500m, na.rm=T)

data[is.na(data$YrsSinceDi),]$YrsSinceDi <- 100
data[is.na(data$ogsi),]$ogsi <- mean(data$ogsi, na.rm=T)
data[is.na(data$Cover_2m_max_all),]$Cover_2m_max_all <- mean(data$Cover_2m_max_all, na.rm=T)
data[is.na(data$Cover_4m_16m),]$Cover_4m_16m <- mean(data$Cover_4m_16m, na.rm=T)
data[is.na(data$Cover_4m_16m_all),]$Cover_4m_16m_all <- mean(data$Cover_4m_16m_all, na.rm=T)
data[is.na(data$p25),]$p25 <- mean(data$p25, na.rm=T)
data[is.na(data$p25_all),]$p25_all <- mean(data$p25_all, na.rm=T)
data[is.na(data$p95),]$p95 <- mean(data$p95, na.rm=T)
data[is.na(data$p95_all),]$p95_all <- mean(data$p95_all, na.rm=T)

#add this?
data$elev2 <- c(scale(data$elev)^2)

####################
#check correlation between covariates
#calculate correlation of variables
covar.num <- dplyr::select_if(data, is.numeric)
correlation <- cor(covar.num)
# write.csv(correlation, file="Data/RSF/covar_correlation_25pointsper.csv")

table(data[,c("AnimlID","used")])

#          used
# AnimlID     0    1
# SG-001   2125   85
# SG-001_2 2300   92
# SG-002    500   20
# SG-003   1200   48
# SG-004    525   21
# SG-005   1450   58
# SG-006   4150  166
# SG-007    750   30
# SG-008   3475  139
# SG-009   5475  219
# SG-010   2500  100
# SG-011   2850  114
# SG-013   1050   42
# SG-014    775   31
# SG-015   1750   70
# SG-017    775   31
# SG-018    775   31
# SG-019   3025  121
# SG-020   1000   40
# SG-021    550   22

####################
#plot rug plots vs elevation

ggplot(data=data, aes(x=elev, y=used)) + geom_jitter(width=0, alpha=0.5) + theme_bw(base_size=20) + facet_wrap(~AnimlID)
# ggsave(filename="Figures/elevation_jitter.tiff", height=10, width=10, units="in", dpi=400, compression="lzw")

####################
#remove animals with not enough information
# data <- data[!data$AnimlID %in% c('SG-023','SG-025'),] #remove animals with <10 locations
data <- data[!data$AnimlID %in% c('SG-012','SG-016','SG-022','SG-023','SG-025'),] #remove animals with <20 locations
# 
# data <- data[data$AnimlID %in% c('SG-006','SG-008','SG-009','SG-010','SG-011','SG-019'),] #animals with >100 locations

###############
#models
###############
#run models for all data together, regardless of animal
#get global and then run for each animal individually

family=binomial(link="logit")

null <- bglmer(used ~ (1|AnimlID), data=data, family=family)


#########
#disturbance
#scaled and centered
#########
d.r <- bglmer(used ~ scale(dist.road) + (1|AnimlID), data=data, family=family)
d.logr <- bglmer(used ~ scale(dist.road) + log(dist.road + 0.01) + (1|AnimlID), data=data, family=family)

d.l <- bglmer(used ~ scale(dist.logging) + (1|AnimlID), data=data, family=family)
d.logl <- bglmer(used ~ scale(dist.logging) + log(dist.logging + 0.01) + (1|AnimlID), data=data, family=family)

y <- bglmer(used ~ scale(YrsSinceDi) + (1|AnimlID), data=data, family=family)
logy <- bglmer(used ~ scale(YrsSinceDi) + log(YrsSinceDi) + (1|AnimlID), data=data, family=family)
o <- bglmer(used ~ scale(ogsi) + (1|AnimlID), data=data, family=family)
logo <- bglmer(used ~ scale(ogsi) + log(ogsi+0.01) + (1|AnimlID), data=data, family=family)

#p_logged
r100 <- bglmer(used ~ scale(percentage_inside.r100) + (1|AnimlID), data=data, family=family)
r500 <- bglmer(used ~ scale(percentage_inside.r500) + (1|AnimlID), data=data, family=family)
r1000 <- bglmer(used ~ scale(percentage_inside.r1000) + (1|AnimlID), data=data, family=family)
r5000 <- bglmer(used ~ scale(percentage_inside.r5000) + (1|AnimlID), data=data, family=family)

#p_>ogsi80 (p_mature and oldgrowth)
# r100.o80 <- bglmer(used ~ scale(percentage_inside.r100.o80), data=data, family=family)
# r500.o80 <- bglmer(used ~ scale(percentage_inside.r500.o80), data=data, family=family)
# r1000.o80 <- bglmer(used ~ scale(percentage_inside.r1000.o80), data=data, family=family)
# r5000.o80 <- bglmer(used ~ scale(percentage_inside.r5000.o80), data=data, family=family)

#p_>ogsi200 (p_oldgrowth)
r100.o200 <- bglmer(used ~ scale(percentage_inside.r100.o200) + (1|AnimlID), data=data, family=family)
r500.o200 <- bglmer(used ~ scale(percentage_inside.r500.o200) + (1|AnimlID), data=data, family=family)
r1000.o200 <- bglmer(used ~ scale(percentage_inside.r1000.o200) + (1|AnimlID), data=data, family=family)
r5000.o200 <- bglmer(used ~ scale(percentage_inside.r5000.o200) + (1|AnimlID), data=data, family=family)

#p_ogsi80<x<ogsi200 (p_mature)
mature.r100 <- bglmer(used ~ scale(p_mature.r100) + (1|AnimlID), data=data, family=family)
mature.r500 <- bglmer(used ~ scale(p_mature.r500) + (1|AnimlID), data=data, family=family)
mature.r1000 <- bglmer(used ~ scale(p_mature.r1000) + (1|AnimlID), data=data, family=family)
mature.r5000 <- bglmer(used ~ scale(p_mature.r5000) + (1|AnimlID), data=data, family=family)

AICctab(null, d.r, d.logr, d.l, d.logl, y, logy, o, logo,
        r100, r500, r1000, r5000, r100.o200, r500.o200, r1000.o200, r5000.o200,
        mature.r100, mature.r500, mature.r1000, mature.r5000, nobs=nrow(data))

#              dAICc df
# d.logr         0.0 4  *
# d.r           19.7 3 
# d.l          128.0 3  correlated w/ dist r by 0.445
# d.logl       129.6 4 
# r1000        140.8 3  *
# r500         149.7 3 
# r100         155.8 3 
# r5000        163.1 3 
# mature.r100  163.3 3  *
# y            166.8 3  correlated with r1000 -0.45
# logy         166.8 4 
# r500.o200    176.8 3  correlated * -0.39 with r1000
# mature.r500  177.8 3 
# r1000.o200   179.8 3 

# null         182.0 2 
# mature.r1000 183.7 3 
# r5000.o200   183.8 3 
# o            184.0 3 
# r100.o200    184.0 3 
# mature.r5000 184.0 3 
# logo         185.1 4

glob.d <- bglmer(used ~ scale(dist.road) + log(dist.road + 0.01) + #scale(dist.logging) +
                   scale(percentage_inside.r1000) + scale(p_mature.r100) + #scale(percentage_inside.r500.o200) + # scale(YrsSinceDi) + 
                   (1|AnimlID), data=data, family=family, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(glob.d)
(coef.d <- data.frame(coef(summary(glob.d))))
ggplot(data=coef.d) + geom_point(aes(x=Estimate, y=row.names(coef.d))) + geom_errorbar(aes(y=row.names(coef.d), xmin=Estimate-Std..Error, xmax=Estimate+Std..Error), width=0.1, lwd=1) + theme_bw(base_size=20)
dotplot(ranef(glob.d, condVar=T)) #plot random effects of individual animals

#with random slope?
# d.logr.r <- bglmer(used ~ scale(dist.road) + log(dist.road + 0.01) + (scale(dist.road)|AnimlID), data=data, family=family)
# d.l.r <- bglmer(used ~ scale(dist.logging) + (scale(dist.logging)|AnimlID), data=data, family=family)
# r1000.r <- bglmer(used ~ scale(percentage_inside.r1000) + (scale(percentage_inside.r1000)|AnimlID), data=data, family=family)
# mature.r100.r <- bglmer(used ~ scale(p_mature.r100) + (scale(p_mature.r100)|AnimlID), data=data, family=family)

AICctab(null, d.logr, d.logr.r, nobs=nrow(data)) #don't use .r
AICctab(null, d.l, d.l.r, nobs=nrow(data))
#      dAICc df
# d.l.r  0.0  5 
# d.l    7.1  3 
# null  61.1  2 
AICctab(null, r1000, r1000.r, nobs=nrow(data))
#         dAICc df
# r1000.r   0.0 5 
# r1000    75.9 3 
# null    117.0 2 
AICctab(null, mature.r100, mature.r100.r, nobs=nrow(data))
#              dAICc df
# mature.r100.r  0.0  5 
# mature.r100    5.7  3 
# null          24.4  2 

dotplot(ranef(d.l.r, condVar=T)) #plot random effects of individual animals
dotplot(ranef(r1000.r, condVar=T)) #plot random effects of individual animals
dotplot(ranef(mature.r100.r, condVar=T)) #plot random effects of individual animals

#########
#predation
#########

p_cancov <- bglmer(used ~ scale(cancov_2017) + (1|AnimlID), data=data, family=family)
p_cancov_con <- bglmer(used ~ scale(cancov_con_2017) + (1|AnimlID), data=data, family=family)
p_cancov_hdw <- bglmer(used ~ scale(cancov_hdw_2017) + (1|AnimlID), data=data, family=family)

p_Cover_4m_16m_all <- bglmer(used ~ scale(Cover_4m_16m_all) + (1|AnimlID), data=data, family=family)
p_Cover_4m_16m <- bglmer(used ~ scale(Cover_4m_16m) + (1|AnimlID), data=data, family=family)
p_Cover_2m_max_all <- bglmer(used ~ scale(Cover_2m_max_all) + (1|AnimlID), data=data, family=family)

p_ht <- bglmer(used ~ scale(canht) + (1|AnimlID), data=data, family=family)
p_stndhgt <- bglmer(used ~ scale(stndhgt_2017) + (1|AnimlID), data=data, family=family)

p_p95_all <- bglmer(used ~ scale(p95_all) + (1|AnimlID), data=data, family=family)
p_p95 <- bglmer(used ~ scale(p95) + (1|AnimlID), data=data, family=family)
p_p25_all <- bglmer(used ~ scale(p25_all) + (1|AnimlID), data=data, family=family)
p_p25 <- bglmer(used ~ scale(p25) + (1|AnimlID), data=data, family=family)

p_tri <- bglmer(used ~ scale(tri) + (1|AnimlID), data=data, family=family)

AICctab(null, p_cancov, p_cancov_con, p_cancov_hdw,
        p_Cover_4m_16m_all, p_Cover_4m_16m, p_Cover_2m_max_all,
        p_ht, p_stndhgt,
        p_p95_all, p_p95, p_p25_all, p_p25, p_tri, nobs=nrow(data))

#                    dAICc df
# p_tri                0.0 3 *
# p_p25               86.8 3 *
# p_p25_all           87.9 3 
# p_Cover_4m_16m_all  95.8 3   correlated with p25 -0.69
# p_ht                98.4 3   correlated with p25 0.51
# p_Cover_4m_16m      99.8 3 
# p_cancov_hdw       116.2 3 

# null               119.3 2 
# p_stndhgt          120.7 3 
# p_cancov           120.8 3 
# p_cancov_con       121.2 3 
# p_p95              121.2 3 
# p_p95_all          121.3 3 
# p_Cover_2m_max_all 121.3 3 

glob.p <- bglmer(used ~ scale(tri) + scale(p25) + (1|AnimlID), data=data, family=family)
summary(glob.p)
coef(summary(glob.p))
dotplot(ranef(glob.p, condVar=T)) #plot random effects of individual animals

#with random slope?
# p_tri.r <- bglmer(used ~ scale(tri) + (scale(tri)|AnimlID), data=data, family=family)
# p_p25.r <- bglmer(used ~ scale(p25) + (scale(p25)|AnimlID), data=data, family=family)
# p_Cover_4m_16m_all.r <- bglmer(used ~ scale(Cover_4m_16m_all) + (scale(Cover_4m_16m_all)|AnimlID), data=data, family=family)
# dotplot(ranef(p_tri.r, condVar=T)) #plot random effects of individual animals
# dotplot(ranef(p_p25.r, condVar=T)) #plot random effects of individual animals
# dotplot(ranef(p_Cover_4m_16m_all.r, condVar=T)) #plot random effects of individual animals

AICctab(null, p_tri, p_tri.r, nobs=nrow(data))
#         dAICc df
# p_tri.r   0.0 5 
# p_tri    25.2 3 
# null    144.5 2 
AICctab(null, p_p25, p_p25.r, nobs=nrow(data)) #don't use random slope for p25
AICctab(null, p_Cover_4m_16m_all, p_Cover_4m_16m_all.r, nobs=nrow(data)) #don't use random slope

#########
#Resource
#########

r_tpi <- bglmer(used ~ scale(tpi) + (1|AnimlID), data=data, family=family)
r_tpi.250m <- bglmer(used ~ scale(tpi.250m) + (1|AnimlID), data=data, family=family)
r_tpi.500m <- bglmer(used ~ scale(tpi.500m) + (1|AnimlID), data=data, family=family)
r_tpi.1000m <- bglmer(used ~ scale(tpi.1000m) + (1|AnimlID), data=data, family=family)

r_d.w <- bglmer(used ~ scale(dist.waterbody) + (1|AnimlID), data=data, family=family)
r_d.logw <- bglmer(used ~ scale(dist.waterbody) + log(dist.waterbody + 0.01) + (1|AnimlID), data=data, family=family)
r_d.w2 <- bglmer(used ~ scale(dist.waterbody) + scale(dist.waterbody^2) + (1|AnimlID), data=data, family=family)

r_d.s <- bglmer(used ~ scale(dist.stream) + (1|AnimlID), data=data, family=family)
# d.s <- bglmer(used ~ dist.stream + (1|AnimlID), data=data.s, family=family)
r_d.logs <- bglmer(used ~ scale(dist.stream) + log(dist.stream + 0.01) + (1|AnimlID), data=data, family=family)

#landsat
r_b2 <- bglmer(used ~ scale(B2_20180717) + (1|AnimlID), data=data, family=family)
r_b3 <- bglmer(used ~ scale(B3_20180717) + (1|AnimlID), data=data, family=family)
r_b4 <- bglmer(used ~ scale(B4_20180717) + (1|AnimlID), data=data, family=family)
r_b5 <- bglmer(used ~ scale(B5_20180717) + (1|AnimlID), data=data, family=family)
r_b6 <- bglmer(used ~ scale(B6_20180717) + (1|AnimlID), data=data, family=family)
r_b7 <- bglmer(used ~ scale(B7_20180717) + (1|AnimlID), data=data, family=family)

r_agedom <- bglmer(used ~ scale(age_dom_2017) + (1|AnimlID), data=data, family=family)

r_sbph <- bglmer(used ~ scale(sbph_ge_25_2017) + (1|AnimlID), data=data, family=family)
r_stph <- bglmer(used ~ scale(stph_ge_25_2017) + (1|AnimlID), data=data, family=family)
r_svph <- bglmer(used ~ scale(svph_ge_25_2017) + (1|AnimlID), data=data, family=family)

r_tph <- bglmer(used ~ scale(tph_ge_3_2017) + (1|AnimlID), data=data, family=family) #this should be in predation
r_tphc <- bglmer(used ~ scale(tphc_ge_3_2017) + (1|AnimlID), data=data, family=family)
r_tphh <- bglmer(used ~ scale(tphh_ge_3_2017) + (1|AnimlID), data=data, family=family)

r_ddi <- bglmer(used ~ scale(ddi_2017) + (1|AnimlID), data=data, family=family)
r_sdi <- bglmer(used ~ scale(sdi_reineke_2017) + (1|AnimlID), data=data, family=family)

AICctab(null,
        r_tpi, r_tpi.250m, r_tpi.500m, r_tpi.1000m,
        r_d.w, r_d.logw, r_d.w2,
        r_d.s, r_d.logs,
        r_b2, r_b3, r_b4, r_b5, r_b6, r_b7, r_agedom,
        r_sbph, r_stph, r_svph,
        r_tph, r_tphc, r_tphh,
        r_ddi, r_sdi,
        nobs=nrow(data))

#             dAICc df
# r_tpi.1000m   0.0 3 *
# r_tpi.500m   34.7 3 
# r_tpi.250m   66.8 3 
# r_d.logw     78.1 4 #dont use this since only blue river reservoir 
# r_d.w2       81.2 4 #dont use this since only blue river reservoir 
# r_d.logs     83.8 4 correlated w/ tpi1000 0.46
# r_d.s        85.1 3 
# r_stph       86.4 3 *
# r_b4         91.9 3 * correlated w/ dist.waterbody -0.386, correlated w/ dist.stream 0.379
# r_b3         92.8 3 
# r_sdi        94.9 3 * correlated w/ stph 0.36688, correlated with b4 by -0.39
# r_agedom     95.2 3 correlated with sdi by 0.5798
# r_sbph       96.2 3 
# r_b5         97.1 3 correlated w/ stph -0.338, correlated with sdi by -0.43776

# null         98.2 2 
# r_svph       98.3 3 
# r_b2         98.5 3 
# r_b7         98.9 3 
# r_tphh       99.8 3 
# r_tphc      100.1 3 
# r_tpi       100.1 3 
# r_d.w       100.1 3 
# r_ddi       100.2 3 
# r_b6        100.2 3 
# r_tph       100.2 3 

glob.r <- bglmer(used ~ scale(tpi.1000m) + 
                   # scale(dist.waterbody) + log(dist.waterbody + 0.01) +
                   scale(dist.stream) + log(dist.stream + 0.01) +
                   scale(stph_ge_25_2017) + scale(B4_20180717) + 
                   scale(sdi_reineke_2017) + 
                   (1|AnimlID), data=data, family=family)
summary(glob.r)
coef(summary(glob.r))
dotplot(ranef(glob.r, condVar=T)) #plot random effects of individual animals

#########
#Thermal tolerance
#########
t_e <- bglmer(used ~ scale(elev) + (1|AnimlID), data=data, family=family)
t_e2 <- bglmer(used ~ scale(elev) + scale(elev^2) + (1|AnimlID), data=data, family=family)

t_abam <- bglmer(used ~ scale(abam_ba_2017) + (1|AnimlID), data=data, family=family)
t_acma <- bglmer(used ~ scale(acma3_ba_2017) + (1|AnimlID), data=data, family=family)
t_psme <- bglmer(used ~ scale(psme_ba_2017) + (1|AnimlID), data=data, family=family)
t_tshe <- bglmer(used ~ scale(tshe_ba_2017) + (1|AnimlID), data=data, family=family)

t_ba <- bglmer(used ~ scale(ba_ge_3_2017) + (1|AnimlID), data=data, family=family)
t_bac <- bglmer(used ~ scale(bac_ge_3_2017) + (1|AnimlID), data=data, family=family)
t_bah <- bglmer(used ~ scale(bah_ge_3_2017) + (1|AnimlID), data=data, family=family)

t_aspect <- bglmer(used ~ scale(aspect) + (1|AnimlID), data=data, family=family)
t_east <- bglmer(used ~ eastness + (1|AnimlID), data=data, family=family)
t_north <- bglmer(used ~ northness + (1|AnimlID), data=data, family=family)


AICctab(null, t_e, t_e2,
        t_abam, t_acma, t_psme, t_tshe,
        t_ba, t_bac, t_bah,
        t_aspect, t_east, t_north,
        nobs=nrow(data))


#         dAICc df
# t_north   0.0  3  *
# t_e2     25.6  4  *
# t_e      33.0  3 
# t_psme   36.4  3  *
# t_acma   40.0  3  *
# t_aspect 42.7  3 
# t_bah    42.8  3  #correlated w/ acma 0.75
# t_east   42.8  3 
# t_tshe   44.2  3  *
# t_ba     45.3  3 
# t_bac    46.9  3 

# null     47.5  2 
# t_abam   49.4  3 

glob.t <- bglmer(used ~ scale(elev) + scale(elev2) + northness + scale(acma3_ba_2017) + scale(psme_ba_2017) + scale(tshe_ba_2017) + (1|AnimlID), data=data, family=family)

summary(glob.t)
coef(summary(glob.t))
dotplot(ranef(glob.t, condVar=T)) #plot random effects of individual animals


###############
#run global models
###############

glob.d <- bglmer(used ~ scale(dist.road) + log(dist.road + 0.01) + #scale(dist.logging) +
                   scale(percentage_inside.r1000) + scale(p_mature.r100) + #scale(percentage_inside.r500.o200) + # scale(YrsSinceDi) + 
                   (1|AnimlID), data=data, family=family, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

glob.p <- bglmer(used ~ scale(tri) + scale(p25) + (1|AnimlID), data=data, family=family) #scale(canht) +

glob.r <- bglmer(used ~ scale(tpi.1000m) + 
                   scale(dist.stream) + log(dist.stream + 0.01) +
                   # scale(dist.waterbody) + log(dist.waterbody + 0.01) +
                   scale(stph_ge_25_2017) + scale(B4_20180717) +
                   scale(sdi_reineke_2017) + 
                   (1|AnimlID), data=data, family=family)

glob.t <- bglmer(used ~ scale(elev) + scale(elev2) + northness + scale(acma3_ba_2017) + scale(psme_ba_2017) + scale(tshe_ba_2017) + (1|AnimlID), data=data, family=family)

AICctab(null, glob.d, glob.p, glob.r, glob.t, logLik=T, weights=T, base=T)

#         logLik  AICc    dLogLik dAICc   df weight
# glob.d -6162.2 12336.3   111.8     0.0 6  1     
# glob.p -6200.7 12409.4    73.2    73.1 4  <0.001
# glob.r -6209.3 12434.5    64.7    98.2 8  <0.001
# glob.t -6226.9 12469.8    47.0   133.5 8  <0.001
# null   -6273.9 12551.9     0.0   215.6 2  <0.001

AICctab(null,
        d.logr, r1000, r500.o200, mature.r100,
        r_tpi.1000m,
        r_d.logs, r_d.s, r_b4, r_b5, r_agedom, r_sdi,
        p_p25, p_tri,
        t_e, t_e2, t_acma, t_psme, t_tshe, t_north,
        nobs=nrow(data))

#             dAICc df
# d.logr        0.0 4 
# p_tri        62.7 3 
# r_tpi.1000m  83.8 3 
# t_north     134.5 3 
# r1000       140.8 3 
# p_p25       149.6 3 
# t_e2        160.1 4 
# mature.r100 163.3 3 
# t_e         167.5 3 
# r_d.logs    167.6 4 #correlated with e 0.4159
# r_d.s       168.8 3 
# t_psme      170.9 3 
# t_acma      174.5 3 
# r_b4        175.7 3 
# r500.o200   176.8 3 
# r_sdi       178.7 3 
# t_tshe      178.7 3 
# r_agedom    178.9 3 
# r_b5        180.8 3 
# null        182.0 2 

#random intercept model, no random slope
#with dist.stream
global <- bglmer(used ~ 
                   scale(dist.road) + log(dist.road + 0.01) + #disturbance
                   scale(percentage_inside.r1000) + scale(p_mature.r100) +
                   
                   scale(tri) + scale(p25) + #predation
                   
                   scale(tpi.1000m) + #resource
                   # scale(dist.waterbody) + log(dist.waterbody + 0.01) +
                   scale(dist.stream) + log(dist.stream + 0.01) +
                   scale(stph_ge_25_2017) + scale(B4_20180717) +
                   #scale(sdi_reineke_2017) + #correlated with psme and tshe
                   
                   scale(elev) + I(scale(elev)^2) + northness + scale(acma3_ba_2017) + scale(psme_ba_2017) + scale(tshe_ba_2017) + #thermal #elev2 is correlated with r1000 by -0.41, elev is correlated w/ tpi.1000m by 0.449; sdi and tshe correlated by 0.577
                   (1|AnimlID), data=data, family=family, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

#no dist.stream
global2 <- bglmer(used ~ 
                    scale(dist.road) + log(dist.road + 0.01) + #disturbance
                    scale(percentage_inside.r1000) + scale(p_mature.r100) +
                    
                    scale(tri) + scale(p25) + #predation
                    
                    scale(tpi.1000m) + #resource
                    # scale(dist.waterbody) + log(dist.waterbody + 0.01) +
                    # scale(dist.stream) + log(dist.stream + 0.01) +
                    scale(stph_ge_25_2017) + scale(B4_20180717) +
                    #scale(sdi_reineke_2017) + #correlated with psme and tshe
                    
                    scale(elev) + I(scale(elev)^2) + northness + scale(acma3_ba_2017) + scale(psme_ba_2017) + scale(tshe_ba_2017) + #thermal #elev2 is correlated with r1000 by -0.41, elev is correlated w/ tpi.1000m by 0.449; sdi and tshe correlated by 0.577
                    (1|AnimlID), data=data, family=family, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


# Random effects:
#   Groups  Name        Variance Std.Dev.
#   AnimlID (Intercept) 0.01318  0.1148  
# Number of obs: 38480, groups:  AnimlID, 20
# 
# Fixed effects:
#                                   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                    -2.808657   0.308645  -9.100  < 2e-16 ***
#   scale(dist.road)               -0.149208   0.055686  -2.679  0.00737 ** 
#   log(dist.road + 0.01)          -0.153028   0.031982  -4.785 1.71e-06 ***
#   scale(percentage_inside.r1000)  0.173448   0.034029   5.097 3.45e-07 ***
#   scale(p_mature.r100)           -0.018526   0.032095  -0.577  0.56378    
#   scale(tri)                     -0.220535   0.031049  -7.103 1.22e-12 ***
#   scale(p25)                     -0.133064   0.030617  -4.346 1.39e-05 ***
#   scale(tpi.1000m)               -0.294881   0.037448  -7.874 3.43e-15 ***
#   scale(dist.stream)             -0.063626   0.052169  -1.220  0.22261    
#   log(dist.stream + 0.01)         0.019510   0.046380   0.421  0.67401    
#   scale(stph_ge_25_2017)         -0.030329   0.032695  -0.928  0.35360    
#   scale(B4_20180717)             -0.115586   0.041294  -2.799  0.00512 ** 
#   scale(elev)                    -0.003931   0.042265  -0.093  0.92590    
#   I(scale(elev)^2)               -0.011032   0.026501  -0.416  0.67719    
#   northness                       0.277860   0.040359   6.885 5.79e-12 ***
#   scale(acma3_ba_2017)           -0.076113   0.032634  -2.332  0.01969 *  
#   scale(psme_ba_2017)             0.018208   0.032185   0.566  0.57158    
#   scale(tshe_ba_2017)             0.025858   0.029755   0.869  0.38483 

#random slope for tri
glob.random1 <- bglmer(used ~ 
                         scale(dist.road) + log(dist.road + 0.01) + #disturbance
                         scale(percentage_inside.r1000) + scale(p_mature.r100) +
                         
                         scale(tri) + scale(p25) + #predation
                         
                         scale(tpi.1000m) + #resource
                         # scale(dist.waterbody) + log(dist.waterbody + 0.01) +
                         scale(dist.stream) + log(dist.stream + 0.01) +
                         scale(stph_ge_25_2017) + scale(B4_20180717) +
                         #scale(sdi_reineke_2017) + #correlated with psme and tshe
                         
                         scale(elev) + I(scale(elev)^2) + northness + scale(acma3_ba_2017) + scale(psme_ba_2017) + scale(tshe_ba_2017) + #thermal #elev2 is correlated with r1000 by -0.41, elev is correlated w/ tpi.1000m by 0.449; sdi and tshe correlated by 0.577
                         (scale(tri)|AnimlID), 
                       data=data, family=family, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


#random slope for p.logged.1km
glob.random2 <- bglmer(used ~ 
                         scale(dist.road) + log(dist.road + 0.01) + #disturbance
                         scale(percentage_inside.r1000) + scale(p_mature.r100) +
                         
                         scale(tri) + scale(p25) + #predation
                         
                         scale(tpi.1000m) + #resource
                         # scale(dist.waterbody) + log(dist.waterbody + 0.01) +
                         scale(dist.stream) + log(dist.stream + 0.01) +
                         scale(stph_ge_25_2017) + scale(B4_20180717) +
                         #scale(sdi_reineke_2017) + #correlated with psme and tshe
                         
                         scale(elev) + I(scale(elev)^2) + northness + scale(acma3_ba_2017) + scale(psme_ba_2017) + scale(tshe_ba_2017) + #thermal #elev2 is correlated with r1000 by -0.41, elev is correlated w/ tpi.1000m by 0.449; sdi and tshe correlated by 0.577
                         (scale(percentage_inside.r1000)|AnimlID), 
                       data=data, family=family, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

#random slope for elev
glob.random10 <- bglmer(used ~ 
                          scale(dist.road) + log(dist.road + 0.01) + #disturbance
                          scale(percentage_inside.r1000) + scale(p_mature.r100) +
                          
                          scale(tri) + scale(p25) + #predation
                          
                          scale(tpi.1000m) + #resource
                          # scale(dist.waterbody) + log(dist.waterbody + 0.01) +
                          scale(dist.stream) + log(dist.stream + 0.01) +
                          scale(stph_ge_25_2017) + scale(B4_20180717) +
                          #scale(sdi_reineke_2017) + #correlated with psme and tshe
                          
                          # scale(elev) + elev2 + 
                          scale(elev) + I(scale(elev)^2) +
                          northness + scale(acma3_ba_2017) + scale(psme_ba_2017) + scale(tshe_ba_2017) + #thermal #elev2 is correlated with r1000 by -0.41, elev is correlated w/ tpi.1000m by 0.449; sdi and tshe correlated by 0.577
                          ((scale(elev) + I(scale(elev)^2))|AnimlID),
                        data=data, family=family, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

#random slope for elev
glob.random10_2 <- bglmer(used ~ 
                            scale(dist.road) + log(dist.road + 0.01) + #disturbance
                            scale(percentage_inside.r1000) + scale(p_mature.r100) +
                            
                            scale(tri) + scale(p25) + #predation
                            
                            scale(tpi.1000m) + #resource
                            # scale(dist.waterbody) + log(dist.waterbody + 0.01) +
                            scale(dist.stream) + log(dist.stream + 0.01) +
                            scale(stph_ge_25_2017) + scale(B4_20180717) +
                            #scale(sdi_reineke_2017) + #correlated with psme and tshe
                            
                            scale(elev) + elev2 +
                            # scale(elev) + I(scale(elev)^2) +
                            northness + scale(acma3_ba_2017) + scale(psme_ba_2017) + scale(tshe_ba_2017) + #thermal #elev2 is correlated with r1000 by -0.41, elev is correlated w/ tpi.1000m by 0.449; sdi and tshe correlated by 0.577
                            ((scale(elev) + elev2)|AnimlID),
                          data=data, family=family, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


summary(glob.random10)
# Random effects:
# Groups  Name             Variance Std.Dev. Corr       
# AnimlID (Intercept)      0.1187   0.3445              
#         scale(elev)      0.3413   0.5842    0.41      
#         I(scale(elev)^2) 0.2237   0.4730   -0.36 -0.19
# Number of obs: 38480, groups:  AnimlID, 20
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                    -2.76042    0.32451  -8.506  < 2e-16 ***
#   scale(dist.road)               -0.13850    0.05631  -2.460  0.01391 *  
#   log(dist.road + 0.01)          -0.15064    0.03212  -4.689 2.74e-06 ***
#   scale(percentage_inside.r1000)  0.16752    0.03657   4.580 4.64e-06 ***
#   scale(p_mature.r100)           -0.01604    0.03269  -0.491  0.62363    
#   scale(tri)                     -0.22432    0.03167  -7.083 1.41e-12 ***
#   scale(p25)                     -0.13349    0.03107  -4.297 1.73e-05 ***
#   scale(tpi.1000m)               -0.29490    0.03909  -7.544 4.57e-14 ***
#   scale(dist.stream)             -0.10258    0.05388  -1.904  0.05694 .  
#   log(dist.stream + 0.01)         0.01721    0.04738   0.363  0.71644    
#   scale(stph_ge_25_2017)         -0.02986    0.03309  -0.902  0.36683    
#   scale(B4_20180717)             -0.09850    0.04223  -2.333  0.01967 *  
#   scale(elev)                     0.02475    0.14968   0.165  0.86868    
#   I(scale(elev)^2)               -0.38811    0.12922  -3.004  0.00267 ** 
#   northness                       0.30005    0.04081   7.352 1.95e-13 ***
#   scale(acma3_ba_2017)           -0.06109    0.03316  -1.842  0.06543 .  
#   scale(psme_ba_2017)             0.01470    0.03246   0.453  0.65072    
#   scale(tshe_ba_2017)             0.02761    0.02982   0.926  0.35452  

#1 tri
#2 percentage_inside.r1000 *
#3 p_mature.r100
#4 dist.road
#5 p25
#6 tpi.1000m
#7 dist.waterbody, dist.stream
#8 stph_ge_25_2017
#9 B4_20180717
#10 elev *
#11 northness
#12 acma
#13 psme
#14 tshe

AICctab(null, glob, glob.random1, glob.random2, glob.random10)

AICctab(null, glob, glob.random, glob.random2, glob.random3, glob.random4, glob.random5, glob.random6,
        glob.random7, glob.random8, glob.random9, glob.random10, glob.random11, glob.random12, glob.random13, glob.random14)
#              dAICc df
# glob.random10   0.0 24
# glob.random2   49.9 21
# glob.random7   62.0 21
# glob.random6  110.4 21
# glob.random13 118.7 21
# glob.random11 119.7 21
# glob.random8  120.5 21
# glob.random1  121.7 21
# glob.random5  126.7 21
# glob.random3  128.1 21

# glob          133.0 19
# glob.random14 135.8 21
# glob.random12 136.8 21
# glob.random9  137.2 21
# glob.random4  229.6 24
# null          583.8 2 

glob <- glob.random10
parameter <- "elevation"

# glob <- glob.random2
# parameter <- "p_logged.r1000"

# glob <- glob
# parameter <- "random_intercept_only"

#fixed effects
summary(glob)
(coef.glob <- data.frame(coef(summary(glob))))
coef.glob$param <- row.names(coef.glob)
coef.glob$param <- gsub(coef.glob$param, pattern="scale\\(", replacement="")
coef.glob$param <- gsub(coef.glob$param, pattern="_ba_2017", replacement="")
coef.glob$param <- gsub(coef.glob$param, pattern="\\)", replacement="")
coef.glob$param <- gsub(coef.glob$param, pattern=" \\+ 0.01", replacement="\\)")
coef.glob$param <- gsub(coef.glob$param, pattern="_ge_25_2017", replacement="")
coef.glob$param <- gsub(coef.glob$param, pattern="_20180717", replacement="")
coef.glob$param <- gsub(coef.glob$param, pattern="percentage_inside", replacement="p_logged")
coef.glob$param <- gsub(coef.glob$param, pattern="acma3", replacement="acma")
coef.glob$param <- gsub(coef.glob$param, pattern="I\\(", replacement="")
coef.glob$param <- gsub(coef.glob$param, pattern="^\\(", replacement="")

# coef.glob$category <- c(NA, rep("disturbance", 4), rep("predation", 2), rep("resource", 5), rep("thermal ", 2), rep("thermal", 4))
coef.glob$category <- c(NA, rep("disturbance", 4), rep("predation", 2), rep("resource", 5), rep("thermal", 6))
coef.glob$param <- toupper(coef.glob$param)
coef.glob$param <- factor(coef.glob$param, levels=coef.glob$param)

ggplot(data=coef.glob[coef.glob$param != toupper("Intercept"),]) + 
  geom_vline(aes(xintercept=0), lty="dashed", col="grey50") + 
  geom_point(aes(x=Estimate, y=param), size=3) + geom_errorbar(aes(y=param, xmin=Estimate-Std..Error, xmax=Estimate+Std..Error), width=0.1, lwd=1) + 
  theme_bw(base_size=20) +theme(panel.grid=element_blank()) + 
  facet_wrap(~category, scales="free", ncol=1) + 
  # force_panelsizes(rows=c(4,2,5,4,2), cols=c(1), respect=F) +
  force_panelsizes(rows=c(4,2,5,6), cols=c(1), respect=F) +
  ylab("Parameter")
ggsave(filename = paste("Figures/rsf/coef_fixed_randomint-slope_", parameter,".tiff", sep=""), height=10, width=8, units="in", dpi=400, compression="lzw")

write.csv(coef.glob, file="Data/RSF/betacoefficients_global_randomslopeelev.csv")

#random effects
dotplot(ranef(glob, condVar=T)) #plot random effects of individual animals

ranef.glob <- as.data.frame(ranef(glob, condVar=T))
ranef.glob$term <- gsub(ranef.glob$term, pattern="^I", replacement="")
ranef.glob$term <- gsub(ranef.glob$term, pattern="scale\\(I", replacement="")
ranef.glob$term <- gsub(ranef.glob$term, pattern="scale", replacement="")
ranef.glob$term <- gsub(ranef.glob$term, pattern="scale\\(percentage_inside", replacement="p_logged")
ranef.glob$term <- gsub(ranef.glob$term, pattern="\\)", replacement="")
ranef.glob$term <- gsub(ranef.glob$term, pattern="\\(", replacement="")

ranef.glob$term <- factor(ranef.glob$term, levels=c("Intercept","elev","elev^2"))

ranef.glob$grp <- factor(ranef.glob$grp, levels=unique(ranef.glob$grp))

ggplot(data=ranef.glob, aes(x=condval, y=grp)) + geom_point(size=3) +
  geom_vline(aes(xintercept=0), lty="dashed", col="grey50") +
  geom_errorbar(aes(xmin=condval-2*condsd, xmax=condval+2*condsd), width=0.1) +
  ylab("") + xlab("") +
  facet_wrap(~term, scales="free_x", nrow=1) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank())
ggsave(filename = paste("Figures/rsf/coef_random_randomint-slope_", parameter,".tiff", sep=""), height=10, width=8, units="in", dpi=400, compression="lzw")

#random effect marginal plot
# glob.pred <- ggpredict(glob, terms=c("elev [all]", "AnimlID [sample=9]"), type="random")
glob.pred <- ggpredict(glob, terms=c("elev [all]", "AnimlID [all]"), type="random")
# plot(glob.pred, ci=F)

# glob.pred <- ggpredict(glob, terms=c("percentage_inside.r1000 [all]", "AnimlID [all]"), type="random")

ggplot(glob.pred, aes(x, predicted, col=group)) + geom_line(lwd=1) + 
  xlab(parameter) + coord_cartesian(ylim=c(0, 0.15)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank())
ggsave(filename=paste("Figures/rsf/global_random_slope_",parameter,".tiff", sep=""),
       height=8, width=10, units="in", dpi=400, compression="lzw")

##########################
#plot marginal plots
##########################
require(dplyr)
require(ggeffects)

plot.marginal <- function(mod, par, par.name)
{
  glob.pred <- ggpredict(mod, terms=c(paste(par, " [all]", sep="")), type="fixed")
  ggplot(glob.pred, aes(x=x, y=predicted)) + geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="grey80") + geom_line(col="black", lwd=1)  +
    xlab(par.name) + ylab("Predicted use") + 
    theme_bw(base_size=20) + theme(panel.grid=element_blank())
  ggsave(filename = paste("Figures/rsf/bglmer/marginal_", par, ".tiff", sep=""), height=6, width=6, units="in", compression="lzw", dpi=400)
}

#########
plot.marginal(mod=glob.random10_2, par="dist.road", par.name="distance to road (m)")
plot.marginal(mod=glob.random10_2, par="percentage_inside.r1000", par.name="percent logged within 1km")
plot.marginal(mod=glob.random10_2, par="p_mature.r100", par.name="percent mature within 0.1km")
plot.marginal(mod=glob.random10_2, par="tpi.1000m", par.name="topographic position index (1km)")
plot.marginal(mod=glob.random10_2, par="tri", par.name="topographic roughness index")
plot.marginal(mod=glob.random10_2, par="p25", par.name="canopy height (m)")
plot.marginal(mod=glob.random10_2, par="dist.stream", par.name="distance to stream (m)")
# plot.marginal(mod=glob.random10, par="dist.waterbody", par.name="distance to waterbody (m)")
plot.marginal(mod=glob.random10_2, par="stph_ge_25_2017", par.name="snag density (trees/ha)")
plot.marginal(mod=glob.random10_2, par="B4_20180717", par.name="LANDSAT B4")
plot.marginal(mod=glob.random10_2, par="northness", par.name="northness")
plot.marginal(mod=glob.random10_2, par="acma3_ba_2017", par.name="Acer macrophyllum basal area")
plot.marginal(mod=glob.random10_2, par="psme_ba_2017", par.name="Pseudotsuga menziesii basal area")
plot.marginal(mod=glob.random10_2, par="tshe_ba_2017", par.name="Tsuga heterophylla basal area")

plot.marginal(mod=glob.random10, par="elev", par.name="elevation (m)")


# mod <- glob.random10_2
# glob.pred <- ggpredict(mod, terms=c("elev [all]", "elev2 [all]"), type="fixed")
mod <- glob.random10
par.name <- "elevation"
glob.pred <- ggpredict(mod, terms=c("elev [all]"), type="fixed")
ggplot(glob.pred, aes(x=x, y=predicted)) + geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="grey80") + geom_line(col="black", lwd=1)  +
  xlab(par.name) + ylab("Predicted probability of occupancy") + 
  theme_bw(base_size=20) + theme(panel.grid=element_blank())
ggsave(filename = paste("Figures/rsf/bglmer/marginal_", par, ".tiff", sep=""), height=6, width=6, units="in", compression="lzw", dpi=400)

summary(global)$coefficients %>% round(digits=2)
#                                Estimate Std. Error z value Pr(>|z|)
# (Intercept)                       -2.80       0.30   -9.24     0.00
# scale(dist.road)                  -0.15       0.06   -2.67     0.01
# log(dist.road + 0.01)             -0.16       0.03   -5.26     0.00
# scale(percentage_inside.r1000)     0.16       0.03    4.92     0.00
# scale(p_mature.r100)              -0.02       0.03   -0.62     0.54
# scale(tri)                        -0.22       0.03   -7.13     0.00
# scale(p25)                        -0.13       0.03   -4.38     0.00
# scale(tpi.1000m)                  -0.30       0.04   -8.17     0.00
# scale(dist.stream)                -0.09       0.05   -1.69     0.09
# log(dist.stream + 0.01)            0.03       0.05    0.55     0.58
# scale(stph_ge_25_2017)            -0.02       0.03   -0.71     0.47
# scale(B4_20180717)                -0.12       0.04   -2.85     0.00
# scale(elev)                        0.02       0.04    0.45     0.65
# I(scale(elev)^2)                  -0.01       0.03   -0.48     0.63
# northness                          0.28       0.04    7.06     0.00
# scale(acma3_ba_2017)              -0.07       0.03   -2.30     0.02
# scale(psme_ba_2017)                0.01       0.03    0.32     0.75
# scale(tshe_ba_2017)                0.02       0.03    0.58     0.56

summary(glob.random10)$coefficients %>% round(digits=2)
#                                Estimate Std. Error z value Pr(>|z|)
# (Intercept)                       -2.77       0.32   -8.75     0.00
# scale(dist.road)                  -0.15       0.06   -2.64     0.01
# log(dist.road + 0.01)             -0.16       0.03   -5.07     0.00
# scale(percentage_inside.r1000)     0.16       0.04    4.49     0.00
# scale(p_mature.r100)              -0.02       0.03   -0.53     0.60
# scale(tri)                        -0.23       0.03   -7.16     0.00
# scale(p25)                        -0.13       0.03   -4.37     0.00
# scale(tpi.1000m)                  -0.30       0.04   -7.86     0.00
# scale(dist.stream)                -0.13       0.05   -2.39     0.02
# log(dist.stream + 0.01)            0.03       0.05    0.56     0.58
# scale(stph_ge_25_2017)            -0.02       0.03   -0.69     0.49
# scale(B4_20180717)                -0.10       0.04   -2.35     0.02
# scale(elev)                        0.08       0.14    0.57     0.57
# I(scale(elev)^2)                  -0.36       0.11   -3.26     0.00
# northness                          0.30       0.04    7.49     0.00
# scale(acma3_ba_2017)              -0.06       0.03   -1.81     0.07
# scale(psme_ba_2017)                0.01       0.03    0.24     0.81
# scale(tshe_ba_2017)                0.02       0.03    0.63     0.53

ranef(glob.random10, condVar=T) #%>% data.frame() %>% mutate_if(is.numeric, round, digits=2)
# $AnimlID
# (Intercept) scale(elev) I(scale(elev)^2)
# SG-001    0.18964747 -0.08436852       0.17671909
# SG-001_2 -0.51258558 -1.09939239       0.01091114
# SG-002    0.05559935  0.50425912       0.19164408
# SG-003    0.04057572 -0.10049668       0.31710679
# SG-004   -0.09233561  0.55654353       0.16195879
# SG-005    0.13158916  0.23909507       0.16674211
# SG-006    0.54062940  0.23076941       0.02526182
# SG-007    0.07337980  0.33598417       0.25095378
# SG-008    0.01174121  0.51884159      -0.88697866
# SG-009   -0.27073317 -0.44865707       0.26768610
# SG-010    0.12748100 -0.17553166      -0.67197843
# SG-011   -0.25996303 -0.35877867       0.42866788
# SG-012   -0.02321762  0.63220952       0.27114500
# SG-013   -0.43366152 -0.21069731       0.63344541
# SG-014    0.30742753  0.09532235      -0.01000918
# SG-015   -0.05633633  0.29452357      -0.56647900
# SG-016    0.24433295  0.61377495      -0.04650704
# SG-017    0.09801541 -0.42918865      -0.28836588
# SG-018    0.47611076 -0.33446934      -0.35739659
# SG-019   -0.06887260  0.81687444      -0.01550121
# SG-020   -0.55619230 -1.00552293       0.27443794
# SG-021    0.01539655 -0.20536520       0.05423944
# SG-022    0.07383486  0.04011063       0.07860697
# SG-023   -0.15463999 -0.38167268       0.45883313
# SG-025    0.14435477  0.09143617      -0.12936883
# 
# with conditional variances for “AnimlID” 