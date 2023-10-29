#############################
#how do home ranges compare to other carnivores?
#why are some home ranges larger than others?
#############################

##############
#load packages
require(plyr)
require(tidyr)

require(ggplot2)
require(ggpubr)

##############
#load data

#skunk info from capture
skunk.info <- read.table("data/AllSkunkCapts.txt", sep=",", header=T)
skunk.info.unique <- ddply(skunk.info, .(Animal.ID, Sex), summarize, 
                           mass=mean(Mass..g., na.rm=T), mass.min=min(Mass..g., na.rm=T), mass.max=max(Mass..g., na.rm=T),
                           length.total=mean(Length.Total..cm., na.rm=T))


#skunk home range info
hr_ca <- read.table("data/HR_ctmm/hr_ca_estimates.txt", sep=",", header=T)

hr <- hr_ca[hr_ca$size == "homerange",]
ca <- hr_ca[hr_ca$size == "corearea",]

#convert from long to wide
hr_ca_wide <- pivot_wider(data=hr_ca, names_from=size, values_from=est, id_cols=c(animal, sex, group))

##############
#1. How do home range sizes compare to other carnivores?

hr_ca[hr_ca$size=="homerange",] %>% group_by(group) %>% dplyr::summarize(min(est), max(est))

hr_ca[hr_ca$size=="homerange",] %>% group_by(group) %>% dplyr::summarize(mean(est), sd(est)/sqrt(length(est)))
hr_ca[hr_ca$size=="corearea",] %>% group_by(group) %>% dplyr::summarize(mean(est), sd(est)/sqrt(length(est)))

otherHRs <- read.csv("data/HR_ctmm/HRsizes_other_carnivores.csv")
otherHRs <- otherHRs[otherHRs$sex != "",]
otherHRs <- otherHRs[otherHRs$year >= 1995,]
otherHRs$HR_km2 <- otherHRs$mean_HR..ha./100

unique(otherHRs$species)

otherHRs$species <- factor(otherHRs$species, levels=c("prairie spotted skunk","western spotted skunk","island spotted skunk","striped skunk",
                                                      "short-tailed weasel","long-tailed weasel","black-footed ferret","Humboldt marten","pacific marten","wolverine",
                                                      "bobcat","mountain lion","coyote","gray wolf"))

ggplot(data=otherHRs, aes(x=mass..kg., y=HR_km2, group=species, col=species, shape=sex)) + 
  geom_errorbar(aes(xmin=mass..kg.-mass_SE, xmax=mass..kg.+mass_SE), width=0) + geom_errorbar(aes(ymin=HR_km2-SE_HR/100, ymax=HR_km2+SE_HR/100), width=0) + geom_point(size=3, stroke=2) + 
  geom_point(data=otherHRs[otherHRs$study == "this study",], aes(x=mass..kg., y=HR_km2, group=species, shape=sex), col="black", size=3, stroke=2) +
  scale_shape_manual(values=c(4,1,2)) + #scale_color_jco() +  
  scale_color_manual(values=c("#A6CEE3","#1F78B4","blue","#0F5257",
                              "#B5CA8D","#8BB174","grey50","#F9B5AC","#EE7674","#6F1A07","#EAC5D8","#6D597A","#FFD151","#C57B57")) +
  # scale_color_brewer(palette="Paired") +
  scale_x_log10(labels=scales::label_comma()) + scale_y_log10(labels=scales::label_comma()) +
  xlab("mass (kg)") + ylab("home range size (km2)") +
  theme_bw(base_size = 20) + theme(panel.grid=element_blank())
# ggsave(filename = "Figures/HR_vs_mass_other_carnivores.tiff", height=8, width=12, units="in", compression="lzw", dpi=400)

##############
#2. Why are some home ranges larger than others?
#extract info from GIS layers using packages raster, rgdal, sf
data <- read.table("data/HR_ctmm/hr_env_covars_2023-05-02.txt", sep=",", header=T)

#convert from raw counts to percent of home range for logging, mature, o80, oldgrowth
data$p.log <- data$logging/data$totalpix*100
data$p.mature <- data$mature/data$totalpix*100
data$p.o80 <- data$o80/data$totalpix*100
data$p.oldgrowth <- data$oldgrowth/data$totalpix*100

######
data[data$animal == "SG-023",] #very very large confidence intervals for HR size so remove from rest of analysis
data <- data[data$animal != "SG-023",]

#################
#linear regressions to determine drivers of home range size
sex <- lm(est ~ sex, data=data)
mass <- lm(est ~ mass, data=data)

sex.mass <- lm(est ~ sex.m + mass + I(sex.m*mass), data=data)

summary(sex)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   11.396      4.625   2.464   0.0216 *
# sexM          14.833      5.781   2.566   0.0173 *

#################
#only female data
# mass.f <- lm(est ~ mass, data=data[data$sex == "F",])
# length.f <- lm(est ~ length.total, data=data[data$sex == "F",])
# logging.f <- lm(homerange ~ p.log, data=data[data$sex == "F",])
# mature.f <- lm(est ~ p.mature, data=data[data$sex == "F",])
# oldgrowth.f <- lm(est ~ p.oldgrowth, data=data[data$sex == "F",])
# emean.f <- lm(est ~ e.mean, data=data[data$sex == "F",])

#################
#only male data
#
locs.m <- lm(est ~ total, data=data[data$sex == "M",])
mass.m <- lm(est ~ mass, data=data[data$sex == "M",])
length.m <- lm(est ~ length.total, data=data[data$sex == "M",])
logging.m <- lm(est ~ p.log, 
                data=data[data$sex == "M",])
mature.m <- lm(est ~ scale(p.mature), 
               data=data[data$sex == "M",])
oldgrowth.m <- lm(est ~ p.oldgrowth, 
                  data=data[data$sex == "M",])
o80.m <- lm(est ~ p.o80, 
            data=data[data$sex == "M",])
emean.m <- lm(est ~ e.mean, data=data[data$sex == "M",])

summary(locs.m)
summary(mass.m)
summary(length.m)
summary(logging.m) #
summary(mature.m) #
summary(oldgrowth.m) #
summary(o80.m) #
summary(emean.m)

require(ggpmisc) #for stat_poly_eq

plot.hr.var <- function(var, varname)
{
  ggplot(data=data[data$sex=="M",], aes(x=eval(parse(text=var)), y=est)) +
    geom_smooth(aes(group=sex, lty=sex), method = "lm", col="black") +
    # geom_text(x=min(data[,var]), y=)
    stat_poly_eq(formula = y ~ scale(x), use_label(c("P")), size=6) +
    stat_poly_eq(formula = y ~ x, use_label(c("eq")), size=6, label.y = 0.85) +
    geom_point(size=3) +
    geom_errorbar(aes(ymin=low, ymax=high), lwd=0.5) +
    xlab(varname) + 
    ylab("home range size (km^2)") +
    theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none")
}

p1 <- plot.hr.var(var="total", varname="number locations")
p2 <- plot.hr.var(var="mass", varname="mass (g)")
p3 <- plot.hr.var(var="length.total", varname="length (cm)")
p4 <- plot.hr.var(var="p.log", varname="percent logged")
p5 <- plot.hr.var(var="p.mature", varname="percent mature")
p6 <- plot.hr.var(var="p.o80", varname="percent OGSI80")
p7 <- plot.hr.var(var="p.oldgrowth", varname="percent old growth")
p8 <- plot.hr.var(var="e.mean", varname="mean elevation (m)")

p0 <- ggarrange(p1, p2, p3, p8, p4, p6, p7, p5)
p0

# ggsave(p0, filename="Figures/homeranges_ctmm_vs_vars_2023-08-19.tiff", height=15, width=15, units="in", compression="lzw", dpi=400)