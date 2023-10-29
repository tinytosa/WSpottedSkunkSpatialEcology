#######################
#ctmm skunk home range overlap
#######################

require(ctmm)

require(ggplot2)
require(ggpubr)

#######################
#load data
a <- read.table("data/HR_ctmm/hr_ca_estimates.txt", sep=",", header=T)

f <- c("SG-005", "SG-007", "SG-008", "SG-009", "SG-015", "SG-016", "SG-017","SG-018","SG-020") #females
nad83z10 <- "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

load("data/HR_ctmm/DATA_all_skunks_error_2022-11-22.rda") #DATA file for all skunks with error ellipses (vhf + gps)
load("data/HR_ctmm/FITS_all_skunks_error_2022-11-16.rda") #FITS file

########################
#calculate KDEs

AKDES <- akde(DATA, FITS, weights=T) #must have verbose=F

######################
#compare male vs. female home ranges
# m1 <- list(AKDES$`SG-002`, AKDES$`SG-006`, AKDES$`SG-011`, AKDES$`SG-013`, AKDES$`SG-021`, AKDES$`SG-023`, AKDES$`SG-025`)
# m2 <- list(AKDES$`SG-001`, AKDES$`SG-001_2`, AKDES$`SG-003`, AKDES$`SG-004`, AKDES$`SG-010`, AKDES$`SG-012`, AKDES$`SG-014`, AKDES$`SG-019`, AKDES$`SG-022`)
# f <- list(AKDES$`SG-005`, AKDES$`SG-007`, AKDES$`SG-008`, AKDES$`SG-009`,AKDES$`SG-015`,AKDES$`SG-016`,AKDES$`SG-017`,AKDES$`SG-018`,AKDES$`SG-020`)

######################
#calculate overlap in KDEs
o <- overlap(AKDES)

o.est <- o$CI
o.est.low <- as.data.frame(o.est[,,1])
o.est.mean <- as.data.frame(o.est[,,2])
o.est.high <- as.data.frame(o.est[,,3])

o.est.low$animal1 <- row.names(o.est.low)
o.est.mean$animal1 <- row.names(o.est.mean)
o.est.high$animal1 <- row.names(o.est.high)

o.est.low <- pivot_longer(as.data.frame(o.est.low), values_to = "overlap.low", names_to = "animal2", cols=starts_with("SG"))
o.est.mean <- pivot_longer(as.data.frame(o.est.mean), values_to = "overlap", names_to = "animal2", cols=starts_with("SG"))
o.est.high <- pivot_longer(as.data.frame(o.est.high), values_to = "overlap.high", names_to = "animal2", cols=starts_with("SG"))

o.est <- merge(o.est.low, o.est.mean, by=c("animal1","animal2"))
o.est <- merge(o.est, o.est.high, by=c("animal1","animal2"))
o.est <- o.est[o.est$animal1 != o.est$animal2,]

#assign sex to animals
o.est$s1 <- "M"
o.est[o.est$animal1 %in% f,]$s1 <- "F"

o.est$s2 <- "M"
o.est[o.est$animal2 %in% f,]$s2 <- "F"

o.est$sexes <- paste(o.est$s1, o.est$s2, sep="") #create combined column of sexes
o.est[o.est$sexes == "FM",]$sexes <- "MF"
table(o.est$sexes)

o.est$animals <- paste(o.est$animal1, o.est$animal2, sep=".") #combine names of animal 1 and animal 2
# write.table(o.est, file="data/HR_ctmm/overlap_estimates_2022-11-20.txt", sep=",")

#############
#plot histograms for each skunk
p.hist <- ggplot(o.est) +
  geom_histogram(aes(x=overlap, fill=sexes), binwidth=0.05) +
  theme_bw(base_size=20) + facet_wrap(~animal1)
p.hist
# ggsave(p.hist, filename="Figures/ctmm_overlap_by_animal.tiff", height=15, width=20, units="in", dpi=330, compression="lzw")

######
#by group instead of by sex
#assign group to animals
o.est <- merge(o.est, a[,c("animal","group")], by.x="animal1", by.y="animal", all.x=T)
o.est <- merge(o.est, a[,c("animal","group")], by.x="animal2", by.y="animal", all.x=T)

o.est$g <- paste(o.est$group.x, o.est$group.y, sep="-")

#plot histograms for each skunk
p.hist <- ggplot(o.est) +
  geom_histogram(aes(x=overlap, fill=g), binwidth=0.05) +
  scale_fill_manual(values=c("grey50","#82204A","#D7B8F3","#D7B8F3","grey50","#558C8C","#D7B8F3","black","grey50")) +
  theme_bw(base_size=20) + facet_wrap(~group.x + animal1) + theme(panel.grid = element_blank())
p.hist
# ggsave(p.hist, filename="Figures/ctmm_overlap_by_animal_bygroup.tiff", height=15, width=20, units="in", dpi=330, compression="lzw")

####################
#find min and max of overlap per animal
o.min <- o.est %>% group_by(animal1, sexes) %>% slice_min(overlap)

o.max <- o.est %>% group_by(animal1, sexes) %>% slice_max(overlap)
# write.table(o.max, file="Data/HR_ctmm/overlap_maximum_2022-11-20.txt", sep=",")

ddply(o.max, .(sexes), summarise, mean(overlap.low), mean(overlap), mean(overlap.high))
#   sexes       ..1       ..2       ..3
# 1    FF 0.6658740 0.7802627 0.8634185
# 2    MF 0.6075320 0.7582945 0.8945954
# 3    MM 0.6453487 0.8293181 0.9419992

#1 is identical distributions, 0 is animals share no area in common

####################
#find overlap by group
ggplot(data=o.est[o.est$g %in% c("F1-F1","M1-M1","M2-M2","F1-M1","F1-M2","M1-M2"),], aes(x=g, y=overlap)) + geom_point(size=3) + theme_bw(base_size=20) + xlab("")
# ggsave(filename="Figures/ctmm_overlap_bygroup.tiff", height=8, width=6, units="in", dpi=400, compression="lzw")

o.est %>% group_by(animal1, g) %>% slice_max(overlap)

o.max <- ddply(o.est, .(animal1, g), summarise, max.low = max(overlap.low), max.est=max(overlap), max.high = max(overlap.high))
ggplot(data=o.max, aes(x=g, y=max.est)) + geom_boxplot() + geom_point(size=3, position=position_dodge2(width=0.25), pch=1) + theme_bw(base_size=20) + xlab("") + ylim(c(0,1))
# ggsave(filename="Figures/ctmm_overlap_bygroup_maxperanimal.tiff", height=8, width=8, units="in", dpi=400, compression="lzw")

ddply(o.max, .(g), summarize, mean(max.est))

summary(aov(data=o.max, max.est ~ g))
TukeyHSD(aov(data=o.max, max.est ~ g))
#                     diff        lwr        upr     p adj
# F1-M1-F1-F1 -0.052487892 -0.2908972 0.18592137 0.9985464
# F1-M2-F1-F1  0.040675056 -0.1977342 0.27908432 0.9997748
# M1-F1-F1-F1 -0.106292285 -0.3611625 0.14857794 0.9161106
# M1-M1-F1-F1 -0.065629489 -0.3204997 0.18924074 0.9956537
# M1-M2-F1-F1  0.014827632 -0.2400426 0.26969786 0.9999999
# M2-F1-F1-F1 -0.038145603 -0.2765549 0.20026366 0.9998609
# M2-M1-F1-F1  0.002243661 -0.2361656 0.24065293 1.0000000
# M2-M2-F1-F1 -0.187708873 -0.4261181 0.05070039 0.2393610
# F1-M2-F1-M1  0.093162948 -0.1452463 0.33157221 0.9409469
# M1-F1-F1-M1 -0.053804393 -0.3086746 0.20106583 0.9989242
# M1-M1-F1-M1 -0.013141597 -0.2680118 0.24172863 1.0000000
# M1-M2-F1-M1  0.067315525 -0.1875547 0.32218575 0.9948346
# M2-F1-F1-M1  0.014342290 -0.2240670 0.25275156 0.9999999
# M2-M1-F1-M1  0.054731553 -0.1836777 0.29314082 0.9980435
# M2-M2-F1-M1 -0.135220981 -0.3736302 0.10318829 0.6697170
# M1-F1-F1-M2 -0.146967341 -0.4018376 0.10790288 0.6498659
# M1-M1-F1-M2 -0.106304545 -0.3611748 0.14856568 0.9160598
# M1-M2-F1-M2 -0.025847424 -0.2807176 0.22902280 0.9999958
# M2-F1-F1-M2 -0.078820659 -0.3172299 0.15958861 0.9778287
# M2-M1-F1-M2 -0.038431395 -0.2768407 0.19997787 0.9998528
# M2-M2-F1-M2 -0.228383929 -0.4667932 0.01002534 0.0707525**
# M1-M1-M1-F1  0.040662796 -0.2296679 0.31099349 0.9999128
# M1-M2-M1-F1  0.121119918 -0.1492108 0.39145062 0.8792434
# M2-F1-M1-F1  0.068146683 -0.1867235 0.32301691 0.9943875
# M2-M1-M1-F1  0.108535946 -0.1463343 0.36340617 0.9064820
# M2-M2-M1-F1 -0.081416587 -0.3362868 0.17345364 0.9820919
# M1-M2-M1-M1  0.080457122 -0.1898736 0.35078782 0.9885930
# M2-F1-M1-M1  0.027483887 -0.2273863 0.28235411 0.9999933
# M2-M1-M1-M1  0.067873150 -0.1869971 0.32274338 0.9945379
# M2-M2-M1-M1 -0.122079383 -0.3769496 0.13279084 0.8343222
# M2-F1-M1-M2 -0.052973235 -0.3078435 0.20189699 0.9990384
# M2-M1-M1-M2 -0.012583971 -0.2674542 0.24228625 1.0000000
# M2-M2-M1-M2 -0.202536505 -0.4574067 0.05233372 0.2287286
# M2-M1-M2-F1  0.040389264 -0.1980200 0.27879853 0.9997864
# M2-M2-M2-F1 -0.149563270 -0.3879725 0.08884600 0.5409695
# M2-M2-M2-M1 -0.189952534 -0.4283618 0.04845673 0.2257632