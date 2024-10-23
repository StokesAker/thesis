library.path <- .libPaths()
library("DHARMa", lib.loc = library.path)

install.packages("glmmTMB")
install.packages("ggpubr")
install.packages("rstatix")
install.packages("broom")
install.packages("ggrepel")
install.packages("pscl")
install.packages("boot")
install.packages("gridExtra")
install.packages("lme4")
install.packages("lsmeans")
install.packages("glmmTMB")
install.packages("TMB", type = "source")
install.packages("DHARMa", dependencies = TRUE, repos = "http://cran.rstudio.com/")
install.packages("GLMMadaptive")
install.packages("vegan")
install.packages("car")
update.packages("car")
install.packages("ggordiplots")
install.packages("permute")
install.packages("lattice")
install.packages("simboot")

library(glmmTMB)
library(ggpubr)
library(rstatix)
library(broom)
library(ggrepel)
library(pscl)
library(boot)
library(gridExtra)
library(lme4)
library(lsmeans)
library(DHARMa)
library(GLMMadaptive)
library(vegan)
library(car)
library(ggordiplots)
library(permute)
library(lattice)
library(simboot)

abund <- read.csv(file = "https://raw.githubusercontent.com/StokesAker/thesis/main/Aker_Abundance_ZerosAdded_2024-2024_v5.csv")

str(abund$site)

abund$time.rep <- as.factor(abund$time.rep)
abund$rep.num <- as.factor(abund$rep.num)
abund$morph <- as.factor(abund$morph)
abund$group <- as.factor(abund$group)

abund2 <- abund %>% 
  dplyr::mutate(treatment = dplyr::recode(treatment,
                                          control = "Control",
                                          nov = "Novaluron",
                                          ace = "Acephate",
                                          bif = "Bifenthrin",
                                          combo = "Combination",
                                          cloth = "Clothianidin")) %>%
  mutate(time.rep = dplyr::recode(time.rep,
                                  "0" = "Pre-treatment",
                                  "1" = "0-7 DAT",
                                  "2" = "7-14 DAT",
                                  "3" = "21-28 DAT",
                                  "4" = "49-56 DAT")) %>%
  mutate(method = dplyr::recode(method, "pitfall" = "Pitfalls", "soilcore" = "Soil Cores"))



abund2$treatment <- factor(abund2$treatment, levels = c("Control","Novaluron","Acephate","Bifenthrin","Clothianidin","Combination"))

str(abund2$site)

###CARABIDAE###
##PiP YEAR/SITE
carabiddt <- filter(abund2, group == "Carabidae", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

carabiddt$tukey <- c("a","ab","a","a","ab",
                     "a","ab","a","a","ab",
                     "a","b","a","a","ab",
                     "a","a","a","a","a",
                     "a","ab","a","a","ab",
                     "a","ab","a","a","b")

carabiddt$tukeytest <- c("1","2","3","4","5",
                         "6","7","8","9","10",
                         "11","12","13","14","15",
                         "16","17","18","19","20",
                         "21","22","23","24","25",
                         "26","27","28","29","30",
                         "31","32","33","34","35",
                         "36","37","38","39","40",
                         "41","42","43","44","45",
                         "46","47","48","49","50",
                         "51","52","53","54","55",
                         "56","57","58","59","60",
                         "61","62","63","64","65",
                         "66","67","68","69","70",
                         "71","72","73","74","75",
                         "76","77","78","79","80",
                         "81","82","83","84","85",
                         "86","87","88","89","90",
                         "91","92","93","94","95",
                         "96","97","98","99","100",
                         "101","102","103","104","105",
                         "106","107","108","109","110",
                         "111","112","113","114","115",
                         "116","117","118","119","120")

carabidplot <- 
  ggplot(carabiddt, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance")+
  theme_minimal()+
  theme(legend.position = c(.20, 1.05),
        legend.justification = c("right","top"),
        legend.box.just = "left",
        legend.margin = margin(2,2,2,2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid(cols = vars(site), rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
carabidplot


##CARABID TT 2023

carabiddtTT23 <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "TT", year ==  "2023") %>%
  aggregate(count ~ treatment + time.rep + rep.num+year, sum) %>%
  group_by(treatment, time.rep, year) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

carabiddtTT23$tukey <- c("a","a","a","a","a",
                     "a","a","a","a","a",
                     "a","a","a","a","a",
                     "a","a","a","a","a",
                     "a","a","a","a","a",
                     "a","a","a","a","a")

carabiddtTT23$tukeytest <- c("1","2","3","4","5",
                         "6","7","8","9","10",
                         "11","12","13","14","15",
                         "16","17","18","19","20",
                         "21","22","23","24","25",
                         "26","27","28","29","30")

carabidplotTT23 <- 
  ggplot(carabiddtTT23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: Toftrees 2023")+
  theme_minimal()+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")

carabidplotTT23

carabiddataTT23 <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "TT", year == "2023")

zerocarabidTT23 <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep), data=carabiddataTT23, family = nbinom2)
summary(zerocarabidTT23)
lsmeans(zerocarabidTT23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocarabidTT23, plot = T)
testZeroInflation(zerocarabidTT23)

###CARABID TT 2024


carabiddtTT24 <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "TT", year ==  "2024") %>%
  aggregate(count ~ treatment + time.rep + rep.num+year, sum) %>%
  group_by(treatment, time.rep, year) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

carabiddtTT24$tukey <- c("a","ac","a","ab","a",
                         "a","a","a","ab","b",
                         "a","a","a","b","ab",
                         "a","b","b","ab","ab",
                         "a","bc","a","ab","ab",
                         "a","bc","a","a","ab")

carabiddtTT24$tukeytest <- c("1","2","3","4","5",
                             "6","7","8","9","10",
                             "11","12","13","14","15",
                             "16","17","18","19","20",
                             "21","22","23","24","25",
                             "26","27","28","29","30")

carabidplotTT24 <- 
  ggplot(carabiddtTT24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: Toftrees 2024")+
  theme_minimal()+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")

carabidplotTT24

carabiddataTT24 <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "TT", year == "2024")

zerocarabidTT24 <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep), data=carabiddataTT24, family = nbinom2)
summary(zerocarabidTT24)
lsmeans(zerocarabidTT24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocarabidTT24, plot = T)
testZeroInflation(zerocarabidTT24)

##CARABID MV 2023
carabiddtMV23 <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "MV", year == "2023") %>%
  aggregate(count ~ treatment + time.rep + rep.num, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

carabiddtMV23$tukey <- c("a","b","a","a","a",
                       "a","b","a","a","a",
                       "a","b","a","a","a",
                       "a","a","a","a","a",
                       "a","b","a","a","a",
                       "a","ab","a","a","a")

carabiddtMV23$tukeytest <- c("1","2","3","4","5",
                           "6","7","8","9","10",
                           "11","12","13","14","15",
                           "16","17","18","19","20",
                           "21","22","23","24","25",
                           "26","27","28","29","30")

carabidplotMV23 <- 
  ggplot(carabiddtMV23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: Mountain View 2023")+
  theme_minimal()+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
carabidplotMV23

carabiddataMV23 <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "MV", year == "2023")

zerocarabidMV23 <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep), data=carabiddataMV23, family = nbinom2)
summary(zerocarabidMV23)
lsmeans(zerocarabidMV23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocarabidMV23, plot = T)
testZeroInflation(zerocarabidMV23)

##CARABID MV 2024
carabiddtMV24 <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "MV", year == "2024") %>%
  aggregate(count ~ treatment + time.rep + rep.num, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

carabiddtMV24$tukey <- c("a","b","a","a","a",
                       "a","b","a","a","a",
                       "a","ab","a","a","a",
                       "a","b","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a")

carabiddtMV24$tukeytest <- c("1","2","3","4","5",
                           "6","7","8","9","10",
                           "11","12","13","14","15",
                           "16","17","18","19","20",
                           "21","22","23","24","25",
                           "26","27","28","29","30")

carabidplotMV24 <- 
  ggplot(carabiddtMV24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: Mountain View 2024")+
  theme_minimal()+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
carabidplotMV24

carabiddataMV24 <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "MV", year == "2024")

zerocarabidMV24 <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep), data=carabiddataMV24, family = nbinom2)
summary(zerocarabidMV24)
lsmeans(zerocarabidMV24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocarabidMV24, plot = T)
testZeroInflation(zerocarabidMV24)

###CARABIDAE TT 2023+2024

carabid2yrdtTT <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "TT") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

carabid2yrdtTT$tukeytest <- c("a","a","a","ac","a",
                            "a","a","ab","a","b",
                            "a","a","a","a","a",
                            "a","a","ab","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","ab",
                            "a","a","a","b","a",
                            "b","a","ab","a","ab",
                            "a","a","a","bc","a",
                            "a","a","ab","a","ab",
                            "a","a","a","bc","a",
                            "a","a","b","a","ab")

carabid2yrplotTT <- 
  ggplot(carabid2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid(rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
carabid2yrplotTT

##BLACK AND WHITE

carabid2yrplotTTbw <- 
  ggplot(carabid2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid(rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))+
  scale_fill_grey(start = 0.1, end = 0.9)
  
carabid2yrplotTTbw

###CARABIDAE MV 2023+2024

carabid2yrdtMV <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "MV") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

carabid2yrdtMV$tukeytest <- c("a","a","ab","b","a",
                              "a","a","a","a","a",
                              "a","a","b","b","a",
                              "a","a","a","a","a",
                              "a","a","ab","ab","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","ab","b","a",
                              "a","a","a","a","a",
                              "a","a","ab","b","a",
                              "a","a","a","a","a")

carabid2yrplotMV <- 
  ggplot(carabid2yrdtMV, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: Mountain View 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
carabid2yrplotMV


###CARABIDAE BOTH YEARS COMBINED SITES SEPARATED

carabiddata2yr <- filter(abund2, group == "Carabidae", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

carabiddata2yr$tukeytest <- c("a","a","b","bc","a",
                                "ab","a","ab","a","a",
                                "a","a","b","bc","a",
                                "ab","a","ab","a","a",
                                "a","a","b","b","a",
                                "ab","a","b","a","a",
                                "a","a","a","a","a",
                                "a","a","ab","a","a",
                                "a","a","b","abc","a",
                                "b","a","ab","a","a",
                                "a","a","b","ac","a",
                                "b","a","a","a","a")

carabidplot2yr <- 
  ggplot(carabiddata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.95, .99),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(rows = vars(site))
carabidplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
carabiddata2yrTT <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "TT")

zerocarabid2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=carabiddata2yrTT, family = nbinom2)
summary(zerocarabid2yrTT)
lsmeans(zerocarabid2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocarabid2yrTT, plot = T)
testZeroInflation(zerocarabid2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
carabiddata2yrMV <- filter(abund2, group == "Carabidae", method == "Pitfalls", site == "MV")

zerocarabid2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=carabiddata2yrMV, family = nbinom2)
summary(zerocarabid2yrMV)
lsmeans(zerocarabid2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocarabid2yrMV, plot = T)
testZeroInflation(zerocarabid2yrMV)


###STAPHYLINIDAE###
staphdt <- filter(abund2, group == "Staphylinidae", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

staphdt$tukey <- c("a","a","a","a","a",
                       "a","a","ab","a","a",
                       "a","a","ab","a","a",
                       "a","a","b","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a")


staphdt$tukeytest <- c("1","2","3","4","5",
                         "6","7","8","9","10",
                         "11","12","13","14","15",
                         "16","17","18","19","20",
                         "21","22","23","24","25",
                         "26","27","28","29","30",
                         "31","32","33","34","35",
                         "36","37","38","39","40",
                         "41","42","43","44","45",
                         "46","47","48","49","50",
                         "51","52","53","54","55",
                         "56","57","58","59","60",
                         "61","62","63","64","65",
                         "66","67","68","69","70",
                         "71","72","73","74","75",
                         "76","77","78","79","80",
                         "81","82","83","84","85",
                         "86","87","88","89","90",
                         "91","92","93","94","95",
                         "96","97","98","99","100",
                         "101","102","103","104","105",
                         "106","107","108","109","110",
                         "111","112","113","114","115",
                         "116","117","118","119","120")


staphplot <- ggplot(staphdt, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Staphylinid Abundance")+
  theme_minimal()+
 theme(legend.position = c(.05, 1.05),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site), rows = vars(year))
staphplot  


##STAPHYLINIDAE TT 2023

staphdtTT23 <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "TT", year == "2023") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep, year) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

staphdtTT23$tukey <- c("a","a","ab","a","a",
                   "a","a","bc","a","a",
                   "a","a","ab","a","a",
                   "a","a","c","a","a",
                   "a","a","ab","a","a",
                   "a","a","a","a","a")

staphdtTT23$tukeytest <- c("1","2","3","4","5",
                       "6","7","8","9","10",
                       "11","12","13","14","15",
                       "16","17","18","19","20",
                       "21","22","23","24","25",
                       "26","27","28","29","30")

staphplotTT23 <- ggplot(staphdtTT23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Staphylinid Abundance: Toftrees 2023")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
staphplotTT23

staphdataTT23 <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "TT", year == "2023")

zerostaphTT23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + rep.num + (1|time.rep), data=staphdataTT23, family = nbinom2)
summary(zerostaphTT23)
lsmeans(zerostaphTT23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerostaphTT23, plot = T)
testZeroInflation(zerostaphTT23)

##STAPHYLINIDAE TT 2024

staphdtTT24 <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "TT", year == "2024") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep, year) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

staphdtTT24$tukey <- c("a","ab","a","ab","a",
                       "a","a","a","ab","a",
                       "a","b","a","b","a",
                       "a","ab","a","a","a",
                       "a","ab","a","ab","a",
                       "a","ab","a","ab","a")

staphdtTT24$tukeytest <- c("1","2","3","4","5",
                           "6","7","8","9","10",
                           "11","12","13","14","15",
                           "16","17","18","19","20",
                           "21","22","23","24","25",
                           "26","27","28","29","30")

staphplotTT24 <- ggplot(staphdtTT24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Staphylinid Abundance: Toftrees 2024")+
  theme_minimal()+
  theme(legend.position = c(.05,1.03),
        legend.justification = c("left","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
staphplotTT24


staphdataTT24 <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "TT", year == "2024")

zerostaphTT24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + rep.num + (1|time.rep), data=staphdataTT24, family = nbinom2)
summary(zerostaphTT24)
lsmeans(zerostaphTT24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerostaphTT24, plot = T)
testZeroInflation(zerostaphTT24)

##STAPHYLINIDAE MV 2023

staphdtMV23 <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "MV", year == "2023") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

staphdtMV23$tukey <- c("a","a","a","a","a",
                     "a","a","a","a","a",
                     "a","a","a","a","a",
                     "a","a","a","a","a",
                     "a","a","a","a","a",
                     "a","a","a","a","a")

staphdtMV23$tukeytest <- c("1","2","3","4","5",
                         "6","7","8","9","10",
                         "11","12","13","14","15",
                         "16","17","18","19","20",
                         "21","22","23","24","25",
                         "26","27","28","29","30")

staphplotMV23 <- ggplot(staphdtMV23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Staphylinid Abundance: Mountain View 2023")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
staphplotMV23

staphdataMV23 <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "MV", year == "2023")

zerostaphMV23 <- glmmTMB(count ~ treatment + time.rep  + rep.num + (1|time.rep), data=staphdataMV23, family = nbinom2)
summary(zerostaphMV23)
lsmeans(zerostaphMV23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerostaphMV23, plot = T)
testZeroInflation(zerostaphMV23)

##STAPHYLINIDAE MV 2024

staphdtMV24 <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "MV", year == "2024") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

staphdtMV24$tukey <- c("b","b","b","b","b",
                       "b","b","b","b","b",
                       "b","b","b","b","b",
                       "a","a","a","a","a",
                       "b","b","b","b","b",
                       "b","b","b","b","b")

staphdtMV24$tukeytest <- c("1","2","3","4","5",
                           "6","7","8","9","10",
                           "11","12","13","14","15",
                           "16","17","18","19","20",
                           "21","22","23","24","25",
                           "26","27","28","29","30")

staphplotMV24 <- ggplot(staphdtMV24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Staphylinid Abundance: Mountain View 2024")+
  theme_minimal()+
  theme(legend.position = c(.85, .95),
        legend.justification = c("left","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
staphplotMV24

staphdataMV24 <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "MV", year == "2024")

zerostaphMV24 <- glmmTMB(count ~ treatment + time.rep +  treatment*time.rep + (1|time.rep), data=staphdataMV24, family = nbinom2)
summary(zerostaphMV24)
lsmeans(zerostaphMV24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerostaphMV24, plot = T)
testZeroInflation(zerostaphMV24)

###STAPHYLINIDAE TT 2023+2024

staph2yrdtTT <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "TT") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

staph2yrdtTT$tukeytest <- c("a","a","a","ab","ab",
                            "a","a","ab","a","a",
                            "a","a","a","b","bc",
                            "a","a","ab","a","a",
                            "a","a","a","ab","ab",
                            "a","a","a","a","a",
                            "a","a","a","ab","c",
                            "a","a","b","a","a",
                            "a","a","a","a","ab",
                            "a","a","ab","a","a",
                            "a","a","a","ab","a",
                            "a","a","ab","a","a")

staph2yrplotTT <- 
  ggplot(staph2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Staphylinid Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
staph2yrplotTT

###STAPHYLINIDAE MV 2023+2024

staph2yrdtMV <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "MV") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

staph2yrdtMV$tukeytest <- c("a","a","a","ab","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","ab","a",
                            "a","a","a","a","a",
                            "a","a","a","b","a",
                            "a","a","a","a","a",
                            "a","a","a","ab","a",
                            "a","a","a","a","a",
                            "a","a","a","ab","a",
                            "a","a","a","a","a")

staph2yrplotMV <- 
  ggplot(staph2yrdtMV, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Staphylinid Abundance: Mountain View 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid(rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
staph2yrplotMV

###sTAPHYLINIDAE BOTH YEARS COMBINED SITES SEPARATED

staphdata2yr <- filter(abund2, group == "Staphylinidae", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

staphdata2yr$tukeytest <- c("a","a","ab","ab","a",
                            "bc","a","a","a","a",
                            "a","a","a","ab","a",
                            "ac","a","a","a","a",
                            "a","a","ab","ab","a",
                            "abc","a","a","a","a",
                            "a","a","b","b","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "b","a","a","a","a",
                            "a","a","ab","ab","a",
                            "b","a","a","a","a")

staphplot2yr <- 
  ggplot(staphdata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Carabid Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.95, 1.10),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
staphplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
staphdata2yrTT <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "TT")

zerostaph2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=staphdata2yrTT, family = nbinom2)
summary(zerostaph2yrTT)
lsmeans(zerostaph2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerostaph2yrTT, plot = T)
testZeroInflation(zerostaph2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
staphdata2yrMV <- filter(abund2, group == "Staphylinidae", method == "Pitfalls", site == "MV")

zerostaph2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=staphdata2yrMV, family = nbinom2)
summary(zerostaph2yrMV)
lsmeans(zerostaph2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerostaph2yrMV, plot = T)
testZeroInflation(zerostaph2yrMV)


###ELATERIDAE###
elatdt <- filter(abund2, group == "Elateridae", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num+year+site, sum) %>%
  group_by(treatment, time.rep, site, year) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

elatdt$tukey <- c("a","abc","a","a","a",
                     "a","abc","a","a","a",
                     "a","a","a","a","a",
                     "a","c","a","a","a",
                     "a","ab","a","a","a",
                     "a","bc","a","a","a")

elatdt$tukeytest <- c("1","2","3","4","5",
                       "6","7","8","9","10",
                       "11","12","13","14","15",
                       "16","17","18","19","20",
                       "21","22","23","24","25",
                       "26","27","28","29","30",
                       "31","32","33","34","35",
                       "36","37","38","39","40",
                       "41","42","43","44","45",
                       "46","47","48","49","50",
                       "51","52","53","54","55",
                       "56","57","58","59","60",
                       "61","62","63","64","65",
                       "66","67","68","69","70",
                       "71","72","73","74","75",
                       "76","77","78","79","80",
                       "81","82","83","84","85",
                       "86","87","88","89","90",
                       "91","92","93","94","95",
                       "96","97","98","99","100",
                       "101","102","103","104","105",
                       "106","107","108","109","110",
                       "111","112","113","114","115",
                       "116","117","118","119","120")

elatplot <- ggplot(elatdt, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Elaterid Abundance")+
  theme_minimal()+
  theme(legend.position = c(.05, 1.05),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site), rows = vars(year))
elatplot

##ELATERIDAE TT 2023

elatdtTT23 <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "TT", year == "2023") %>%
  aggregate(count ~ treatment + time.rep + rep.num, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

elatdtTT23$tukey <- c("a","ac","ab","a","a",
                  "a","abc","ab","a","a",
                  "a","bc","b","a","a",
                  "a","ac","a","a","a",
                  "a","c","b","a","a",
                  "a","ab","ab","a","a")

elatdtTT23$tukeytest <- c("1","2","3","4","5",
                      "6","7","8","9","10",
                      "11","12","13","14","15",
                      "16","17","18","19","20",
                      "21","22","23","24","25",
                      "26","27","28","29","30")

elatplotTT23 <- ggplot(elatdtTT23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Elaterid Abundance: Toftrees 2023")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
elatplotTT23

elatdataTT23 <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "TT", year == "2023")

zeroelatTT23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=elatdataTT23, family = nbinom2)
summary(zeroelatTT23)
lsmeans(zeroelatTT23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroelatTT23, plot = T)
testZeroInflation(zeroelatTT23)


##ELATERIDAE TT 2024

elatdtTT24 <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "TT", year == "2024") %>%
  aggregate(count ~ treatment + time.rep + rep.num, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

elatdtTT24$tukey <- c("a","ab","a","a","a",
                      "a","ab","a","a","a",
                      "a","b","a","a","a",
                      "a","ab","a","a","a",
                      "a","ab","a","a","a",
                      "a","a","a","a","a")

elatdtTT24$tukeytest <- c("1","2","3","4","5",
                          "6","7","8","9","10",
                          "11","12","13","14","15",
                          "16","17","18","19","20",
                          "21","22","23","24","25",
                          "26","27","28","29","30")

elatplotTT24 <- ggplot(elatdtTT24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Elaterid Abundance: Toftrees 2024")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
elatplotTT24

elatdataTT24 <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "TT", year == "2024")

zeroelatTT24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=elatdataTT24, family = nbinom2)
summary(zeroelatTT24)
lsmeans(zeroelatTT24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroelatTT24, plot = T)
testZeroInflation(zeroelatTT24)

##ELATERIDAE MV 2023

elatdtMV23 <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "MV", year == "2023") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

elatdtMV23$tukey <- c("a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a")

elatdtMV23$tukeytest <- c("1","2","3","4","5",
                        "6","7","8","9","10",
                        "11","12","13","14","15",
                        "16","17","18","19","20",
                        "21","22","23","24","25",
                        "26","27","28","29","30")

elatplotMV23 <- ggplot(elatdtMV23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Elaterid Abundance: Mountain View 2023")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
elatplotMV23

elatdataMV23 <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "MV", year == "2023")

zeroelatMV23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=elatdataMV23, family = nbinom2)
summary(zeroelatMV23)
lsmeans(zeroelatMV23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroelatMV23, plot = T)
testZeroInflation(zeroelatMV23)

##ELATERIDAE MV 2024

elatdtMV24 <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "MV", year == "2024") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

elatdtMV24$tukey <- c("b","b","b","a","a",
                      "ab","ab","b","ab","ab",
                      "ab","ab","b","ab","ab",
                      "a","a","a","a","b",
                      "b","b","b","b","ab",
                      "ab","ab","ab","b","a")

elatdtMV24$tukeytest <- c("1","2","3","4","5",
                          "6","7","8","9","10",
                          "11","12","13","14","15",
                          "16","17","18","19","20",
                          "21","22","23","24","25",
                          "26","27","28","29","30")

elatplotMV24 <- ggplot(elatdtMV24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Elaterid Abundance: Mountain View 2024")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
elatplotMV24

elatdataMV24 <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "MV", year == "2024")

zeroelatMV24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=elatdataMV24, family = nbinom2)
summary(zeroelatMV24)
lsmeans(zeroelatMV24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroelatMV24, plot = T)
testZeroInflation(zeroelatMV24)

###ELATERIDAE TT 2023+2024

elat2yrdtTT <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "TT") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

elat2yrdtTT$tukeytest <- c("a","a","ab","ab","ab",
                           "a","a","a","a","a",
                           "a","a","ab","ab","ab",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","b","ab","b",
                           "a","a","a","a","a",
                           "a","a","ab","ab","a",
                           "a","a","a","a","a",
                           "a","a","ab","b","ab",
                           "a","a","a","a","a")

elat2yrplotTT <- 
  ggplot(elat2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Elaterid Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid(rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
elat2yrplotTT

###ELATERIDAE MV 2023+2024

elat2yrdtMV <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "MV") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

elat2yrdtMV$tukeytest <- c("a","_","a","_","a",
                           "_","a","_","a","_",
                           "a","_","a","_","a",
                           "_","a","_","a","_",
                           "a","_","a","_","a",
                           "_","a","_","a","_",
                           "a","_","a","_","a",
                           "_","a","_","a","_",
                           "a","_","a","_","a",
                           "_","a","_","a","_",
                           "a","_","a","_","a",
                           "_","a","_","a","_")

elat2yrplotMV <- 
  ggplot(elat2yrdtMV, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Elaterid Abundance: Mountain View 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
elat2yrplotMV

###ELATERIDAE BOTH YEARS COMBINED SITES SEPARATED

elatdata2yr <- filter(abund2, group == "Elateridae", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

elatdata2yr$tukeytest <- c("a","a","ab","bc","a",
                           "ab","a","a","a","a",
                           "a","a","ab","abc","a",
                           "ab","a","a","a","a",
                           "a","a","a","a","a",
                           "b","a","a","a","a",
                           "a","a","b","b","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "b","a","a","a","a",
                           "a","a","b","b","a",
                           "ab","a","a","a","a")

elatplot2yr <- 
  ggplot(elatdata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Elaterid Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.25, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
elatplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
elatdata2yrTT <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "TT")

zeroelat2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=elatdata2yrTT, family = nbinom2)
summary(zeroelat2yrTT)
lsmeans(zeroelat2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroelat2yrTT, plot = T)
testZeroInflation(zeroelat2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
elatdata2yrMV <- filter(abund2, group == "Elateridae", method == "Pitfalls", site == "MV")

zeroelat2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=elatdata2yrMV, family = nbinom1)
summary(zeroelat2yrMV)
lsmeans(zeroelat2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroelat2yrMV, plot = T)
testZeroInflation(zeroelat2yrMV)

###FORMICIDAE###
formdt <- filter(abund2, group == "Formicidae", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num+year+site, sum) %>%
  group_by(treatment, time.rep, site, year) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

formdt$tukey <- c("a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a")

formdt$tukeytest <- c("1","2","3","4","5",
                      "6","7","8","9","10",
                      "11","12","13","14","15",
                      "16","17","18","19","20",
                      "21","22","23","24","25",
                      "26","27","28","29","30",
                      "31","32","33","34","35",
                      "36","37","38","39","40",
                      "41","42","43","44","45",
                      "46","47","48","49","50",
                      "51","52","53","54","55",
                      "56","57","58","59","60",
                      "61","62","63","64","65",
                      "66","67","68","69","70",
                      "71","72","73","74","75",
                      "76","77","78","79","80",
                      "81","82","83","84","85",
                      "86","87","88","89","90",
                      "91","92","93","94","95",
                      "96","97","98","99","100",
                      "101","102","103","104","105",
                      "106","107","108","109","110",
                      "111","112","113","114","115",
                      "116","117","118","119","120")

formplot <- ggplot(formdt, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Formicid Abundance")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site), rows = vars(year))
formplot

##FORMICIDAE TT 2023

formdtTT23 <- filter(abund2, group == "Formicidae", site == "TT", year == "2023", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

formdtTT23$tukey <- c("a","a","a","a","a",
                  "a","a","a","a","a",
                  "a","a","a","a","a",
                  "a","a","a","a","a",
                  "a","a","a","a","a",
                  "a","a","a","a","a")

formdtTT23$tukeytest <- c("1","2","3","4","5",
                      "6","7","8","9","10",
                      "11","12","13","14","15",
                      "16","17","18","19","20",
                      "21","22","23","24","25",
                      "26","27","28","29","30")

formplotTT23 <- ggplot(formdtTT23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Formicid Abundance: Toftrees 2023")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
formplotTT23


formdataTT23 <- filter(abund2, group == "Formicidae", site == "TT", year == "2023", method == "Pitfalls")

zeroformTT23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=formdataTT23, family = nbinom2)
summary(zeroformTT23)
lsmeans(zeroformTT23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroformTT23, plot = T)
testZeroInflation(zeroformTT23)

##FORMICIDAE TT 2024

formdtTT24 <- filter(abund2, group == "Formicidae", site == "TT", year == "2024", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

formdtTT24$tukey <- c("a","a","a","a","a",
                      "a","a","a","a","a",
                      "a","a","a","a","a",
                      "a","a","a","a","a",
                      "a","a","a","a","a",
                      "a","a","a","a","a")

formdtTT24$tukeytest <- c("1","2","3","4","5",
                          "6","7","8","9","10",
                          "11","12","13","14","15",
                          "16","17","18","19","20",
                          "21","22","23","24","25",
                          "26","27","28","29","30")

formplotTT24 <- ggplot(formdtTT24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Formicid Abundance:Toftrees 2024")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
formplotTT24

formdataTT24 <- filter(abund2, group == "Formicidae", site == "TT", year == "2024", method == "Pitfalls")

zeroformTT24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=formdataTT23, family = nbinom2)
summary(zeroformTT24)
lsmeans(zeroformTT24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroformTT24, plot = T)
testZeroInflation(zeroformTT24)

##FORMICIDAE MV 2023

formdtMV23 <- filter(abund2, group == "Formicidae", site == "MV", year == "2023", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

formdtMV23$tukey <- c("a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a",
                    "a","a","a","a","a")

formdtMV23$tukeytest <- c("1","2","3","4","5",
                        "6","7","8","9","10",
                        "11","12","13","14","15",
                        "16","17","18","19","20",
                        "21","22","23","24","25",
                        "26","27","28","29","30")

formplotMV23 <- ggplot(formdtMV23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Formicid Abundance: Mountain View 2023")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
formplotMV23

formdataMV23 <- filter(abund2, group == "Formicidae", site == "MV", year == "2023", method == "Pitfalls")

zeroformMV23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=formdataMV23, family = nbinom2)
summary(zeroformMV23)
lsmeans(zeroformMV23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroformMV23, plot = T)
testZeroInflation(zeroformMV23)

##FORMICIDAE MV 2024

formdtMV24 <- filter(abund2, group == "Formicidae", site == "MV", year == "2024", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

formdtMV24$tukey <- c("a","a","a","a","a",
                      "a","a","a","a","a",
                      "a","a","a","a","a",
                      "a","a","a","a","a",
                      "a","a","a","a","a",
                      "a","a","a","a","a")

formdtMV24$tukeytest <- c("1","2","3","4","5",
                          "6","7","8","9","10",
                          "11","12","13","14","15",
                          "16","17","18","19","20",
                          "21","22","23","24","25",
                          "26","27","28","29","30")

formplotMV24 <- ggplot(formdtMV24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Formicid Abundance: Mountain View 2024")+
  theme_minimal()+
  theme(legend.position = c(.05, .95),
        legend.justification = c("left","top"),
        legend.box.just = "left",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
formplotMV24

formdataMV24 <- filter(abund2, group == "Formicidae", site == "MV", year == "2024", method == "Pitfalls")

zeroformMV24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=formdataMV24, family = nbinom2)
summary(zeroformMV24)
lsmeans(zeroformMV24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroformMV24, plot = T)
testZeroInflation(zeroformMV24)

###FORMICIDAE TT 2023+2024

form2yrdtTT <- filter(abund2, group == "Formicidae", method == "Pitfalls", site == "TT") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

form2yrdtTT$tukeytest <- c("a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a")

form2yrplotTT <- 
  ggplot(form2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Formicid Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
form2yrplotTT

###FORMICIDAE MV 2023+2024

form2yrdtMV <- filter(abund2, group == "Formicidae", method == "Pitfalls", site == "MV") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

form2yrdtMV$tukeytest <- c("a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a")

form2yrplotMV <- 
  ggplot(form2yrdtMV, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Formicid Abundance: Mountain View 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid(rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
form2yrplotMV


###FORMICIDAE BOTH YEARS COMBINED SITES SEPARATED

formdata2yr <- filter(abund2, group == "Formicidae", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

formdata2yr$tukeytest <- c("a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a")

formplot2yr <- 
  ggplot(formdata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Formicid Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.25, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
formplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
formdata2yrTT <- filter(abund2, group == "Formicidae", method == "Pitfalls", site == "TT")

zeroform2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=formdata2yrTT, family = nbinom2)
summary(zeroform2yrTT)
lsmeans(zeroform2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroform2yrTT, plot = T)
testZeroInflation(zeroform2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
formdata2yrMV <- filter(abund2, group == "Formicidae", method == "Pitfalls", site == "MV")

zeroform2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=formdata2yrMV, family = nbinom2)
summary(zeroform2yrMV)
lsmeans(zeroform2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroform2yrMV, plot = T)
testZeroInflation(zeroform2yrMV)

##COLLEMBOLA###

collemboladt <-filter(abund2, group == "Collembola", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site + year, sum) %>%
  group_by(treatment, time.rep, site , year) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

collemboladt$tukey <- c("a","a","a","a","a",
                        "a","a","a","a","a",
                        "a","a","a","a","a",
                        "a","a","a","a","a",
                        "a","a","a","a","a",
                        "a","a","a","a","a")

collemboladt$tukeytest <- c("1","2","3","4","5",
                      "6","7","8","9","10",
                      "11","12","13","14","15",
                      "16","17","18","19","20",
                      "21","22","23","24","25",
                      "26","27","28","29","30",
                      "31","32","33","34","35",
                      "36","37","38","39","40",
                      "41","42","43","44","45",
                      "46","47","48","49","50",
                      "51","52","53","54","55",
                      "56","57","58","59","60",
                      "61","62","63","64","65",
                      "66","67","68","69","70",
                      "71","72","73","74","75",
                      "76","77","78","79","80",
                      "81","82","83","84","85",
                      "86","87","88","89","90",
                      "91","92","93","94","95",
                      "96","97","98","99","100",
                      "101","102","103","104","105",
                      "106","107","108","109","110",
                      "111","112","113","114","115",
                      "116","117","118","119","120")

collembolaplot <-ggplot(collemboladt, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Collembola Abundance (Pitfalls Only)")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site), rows = vars(year))
collembolaplot

## COLLEMBOLA TT 2023

collemboladtTT23 <-filter(abund2, group == "Collembola", site == "TT", year == "2023", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

collemboladtTT23$tukey <- c("a","a","a","a","a",
                          "a","a","a","a","a",
                          "a","a","a","a","a",
                          "a","a","a","a","a",
                          "a","a","a","a","a",
                          "a","a","a","a","a")

collembolaplotTT23 <-ggplot(collemboladtTT23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Collembola Abundance: Toftrees 2023")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
collembolaplotTT23

collemboladataTT23 <- filter(abund2, group == "Collembola", site == "TT", year == "2023", method == "Pitfalls")

zerocollembolaTT23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=collemboladataTT23, family = nbinom2)
summary(zerocollembolaTT23)
lsmeans(zerocollembolaTT23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocollembolaTT23, plot = T)
testZeroInflation(zerocollembolaTT23)

## COLLEMBOLA TT 2024

collemboladtTT24 <-filter(abund2, group == "Collembola", site == "TT", year == "2024", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

collemboladtTT24$tukey <- c("a","a","a","ab","a",
                            "a","a","a","b","b",
                            "a","a","a","ab","a",
                            "a","a","a","a","ab",
                            "a","a","a","ab","ab",
                            "a","a","a","a","a")

collembolaplotTT24 <-ggplot(collemboladtTT24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Collembola Abundance: Toftrees 2024")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
collembolaplotTT24

collemboladataTT24 <- filter(abund2, group == "Collembola", site == "TT", year == "2024", method == "Pitfalls")

zerocollembolaTT24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=collemboladataTT24, family = nbinom2)
summary(zerocollembolaTT24)
lsmeans(zerocollembolaTT24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocollembolaTT24, plot = T)
testZeroInflation(zerocollembolaTT24)

## COLLEMBOLA MV 2023

collemboladtMV23 <-filter(abund2, group == "Collembola", site == "MV", year == "2023", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

collemboladtMV23$tukey <- c("a","a","a","a","a",
                          "a","a","a","a","a",
                          "a","a","a","a","a",
                          "a","a","a","a","a",
                          "a","a","a","a","a",
                          "a","a","a","a","a")

collembolaplotMV23 <-ggplot(collemboladtMV23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Collembola Abundance: Mountain View 2023")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
collembolaplotMV23

collemboladataMV23 <- filter(abund2, group == "Collembola", site == "MV", year == "2023", method == "Pitfalls")

zerocollembolaMV23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=collemboladataMV23, family = nbinom2)
summary(zerocollembolaMV23)
lsmeans(zerocollembolaMV23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocollembolaMV23, plot = T)
testZeroInflation(zerocollembolaMV23)

## COLLEMBOLA MV 2024

collemboladtMV24 <-filter(abund2, group == "Collembola", site == "MV", year == "2024", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

collemboladtMV24$tukey <- c("a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a")


collemboladtMV24$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

collembolaplotMV24 <-ggplot(collemboladtMV24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Collembola Abundance: Mountain View 2024")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
collembolaplotMV24

collemboladataMV24 <- filter(abund2, group == "Collembola", site == "MV", year == "2024", method == "Pitfalls")

zerocollembolaMV24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=collemboladataMV23, family = nbinom2)
summary(zerocollembolaMV24)
lsmeans(zerocollembolaMV24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocollembolaMV24, plot = T)
testZeroInflation(zerocollembolaMV24)

###COLLEMBOLA TT 2023+2024

coll2yrdtTT <- filter(abund2, group == "Collembola", method == "Pitfalls", site == "TT") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

coll2yrdtTT$tukeytest <- c("a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a")
coll2yrplotTT <- 
  ggplot(coll2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Collembola Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
coll2yrplotTT

###COLLEMBOLA MV 2023+2024

coll2yrdtMV <- filter(abund2, group == "Collembola", method == "Pitfalls", site == "MV") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

coll2yrdtMV$tukeytest <- c("a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a",
                           "a","a","a","a","a")
coll2yrplotMV <- 
  ggplot(coll2yrdtMV, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Collembola Abundance: Mountain View 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
coll2yrplotMV

###COLLEMBOLA BOTH YEARS COMBINED SITES SEPARATED

collemboladata2yr <- filter(abund2, group == "Collembola", method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

collemboladata2yr$tukeytest <- c("a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a",
                                 "a","a","a","a","a")
collembolaplot2yr <- 
  ggplot(collemboladata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Collembola Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.25, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
collembolaplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
collemboladata2yrTT <- filter(abund2, group == "Collembola", method == "Pitfalls", site == "TT")

zerocollembola2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=collemboladata2yrTT, family = nbinom2)
summary(zerocollembola2yrTT)
lsmeans(zerocollembola2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocollembola2yrTT, plot = T)
testZeroInflation(zerocollembola2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
collemboladata2yrMV <- filter(abund2, group == "Collembola", method == "Pitfalls", site == "MV")

zerocollembola2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=collemboladata2yrMV, family = nbinom2)
summary(zerocollembola2yrMV)
lsmeans(zerocollembola2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerocollembola2yrMV, plot = T)
testZeroInflation(zerocollembola2yrMV)


###SOIL CORE ACARI###
acaridtSC <-  filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", method == "Soil Cores") %>% 
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year , site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

acaridtSC$tukey <- c(    "a","a","a","a","a",
                         "a","a","a","a","a",
                         "a","a","a","a","a",
                         "a","a","a","a","a",
                         "a","a","a","a","a",
                         "a","a","a","a","a")

acaridtSC$tukeytest <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30",
                            "31","32","33","34","35",
                            "36","37","38","39","40",
                            "41","42","43","44","45",
                            "46","47","48","49","50",
                            "51","52","53","54","55",
                            "56","57","58","59","60",
                            "61","62","63","64","65",
                            "66","67","68","69","70",
                            "71","72","73","74","75",
                            "76","77","78","79","80",
                            "81","82","83","84","85",
                            "86","87","88","89","90",
                            "91","92","93","94","95",
                            "96","97","98","99","100",
                            "101","102","103","104","105",
                            "106","107","108","109","110",
                            "111","112","113","114","115",
                            "116","117","118","119","120")

acariplotSC <-ggplot(acaridtSC, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Soil Cores")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site), rows = vars(year))
acariplotSC

## ACARI TT 2023

acaridtTT23 <-filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", year == "2023", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

acaridtTT23$tukey <-      c("a","ab","ab","a","a",
                            "a","b","ab","a","a",
                            "a","a","a","a","a",
                            "a","ab","b","a","a",
                            "a","ab","ab","a","a",
                            "a","ab","ab","a","a")

acaridtTT23$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

acariplotTT23 <-ggplot(acaridtTT23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Toftrees 2023")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
acariplotTT23

acaridataTT23 <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", year == "2023", method == "Soil Cores")

zeroacariTT23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=acaridataTT23, family = nbinom2)
summary(zeroacariTT23)
lsmeans(zeroacariTT23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroacariTT23, plot = T)
testZeroInflation(zeroacariTT23)

## ACARI TT 2024

acaridtTT24 <-filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", year == "2024", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

acaridtTT24$tukey <- c("a","a","b","a","a",
                       "a","a","ab","a","a",
                       "a","a","ab","a","a",
                       "a","a","ab","a","a",
                       "a","a","b","a","a",
                       "a","a","a","a","a")

acaridtTT24$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

acariplotTT24 <-ggplot(acaridtTT24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Toftrees 2024")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
acariplotTT24

acaridataTT24 <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", year == "2024", method == "Soil Cores")

zeroacariTT24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=acaridataTT24, family = nbinom2)
summary(zeroacariTT24)
lsmeans(zeroacariTT24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroacariTT24, plot = T)
testZeroInflation(zeroacariTT24)

## ACARI MV 2023

acaridtMV23 <-filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", year == "2023", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

acaridtMV23$tukey <- c("a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a")

acaridtMV23$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

acariplotMV23 <-ggplot(acaridtMV23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Mountain View 2023")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
acariplotMV23

acaridataMV23 <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", year == "2023", method == "Soil Cores")

zeroacariMV23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=acaridataMV23, family = nbinom2)
summary(zeroacariMV23)
lsmeans(zeroacariMV23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroacariMV23, plot = T)
testZeroInflation(zeroacariMV23)

## ACARI MV 2024

acaridtMV24 <-filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", year == "2024", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

acaridtMV24$tukey <- c("a","a","a","a","a",
                       "a","a","ab","a","a",
                       "a","a","ab","a","a",
                       "a","b","b","a","a",
                       "a","a","ab","a","a",
                       "a","a","ab","a","a")

acaridtMV24$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

acariplotMV24 <-ggplot(acaridtMV24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Mountain View 2024")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
acariplotMV24

acaridataMV24 <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", year == "2024", method == "Soil Cores")

zeroacariMV24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=acaridataMV24, family = nbinom2)
summary(zeroacariMV24)
lsmeans(zeroacariMV24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroacariMV24, plot = T)
testZeroInflation(zeroacariMV24)

###ACARI TT 2023+2024

acari2yrdtTT <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

acari2yrdtTT$tukeytest <- c("a","a","ab","a","ab",
                            "a","a","a","a","a",
                            "a","a","a","a","ab",
                            "ab","a","a","a","a",
                            "a","a","b","a","b",
                            "ab","a","a","a","a",
                            "a","a","ab","a","a",
                            "ab","a","a","a","a",
                            "a","a","ab","a","ab",
                            "a","a","a","a","a",
                            "a","a","ab","a","ab",
                            "b","a","a","a","a")
acari2yrplotTT <- 
  ggplot(acari2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Acari Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
acari2yrplotTT

###ACARI MV 2023+2024

acari2yrdtMV <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

acari2yrdtMV$tukeytest <- c("a","a","a","b","a",
                            "b","a","a","a","a",
                            "a","a","a","b","a",
                            "ab","a","a","a","a",
                            "a","a","a","b","a",
                            "ab","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","b","a",
                            "ab","a","a","a","a",
                            "a","a","a","b","a",
                            "ab","a","a","a","a")
acari2yrplotMV <- 
  ggplot(acari2yrdtMV, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Acari Abundance: Mountain View 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
acari2yrplotMV

###ACARI BOTH YEARS COMBINED SITES SEPARATED

acaridata2yr <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

acaridata2yr$tukeytest <- c("a","a","b","a","a",
                            "a","a","a","ab","ab",
                            "a","a","ab","a","a",
                            "a","a","a","b","b",
                            "a","a","ab","a","a",
                            "a","a","a","ab","ab",
                            "a","a","a","a","a",
                            "a","a","a","ab","ab",
                            "a","a","ab","a","a",
                            "a","a","a","ab","ab",
                            "a","a","ab","a","a",
                            "a","a","a","a","a")
acariplot2yr <- 
  ggplot(acaridata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Acari Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.25, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
acariplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
acaridata2yrTT <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", method == "Soil Cores", site == "TT")

zeroacari2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=acaridata2yrTT, family = nbinom2)
summary(zeroacari2yrTT)
lsmeans(zeroacari2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroacari2yrTT, plot = T)
testZeroInflation(zeroacari2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
acaridata2yrMV <- filter(abund2, group == "Oribatid Mites" | group == "Non-oribatid Mites" | group == "Predatory Mites", method == "Soil Cores", site == "MV")

zeroacari2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=acaridata2yrMV, family = nbinom2)
summary(zeroacari2yrMV)
lsmeans(zeroacari2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroacari2yrMV, plot = T)
testZeroInflation(zeroacari2yrMV)

###NON-ORIBATID MITES TT 2023
nonoridtTT23 <-filter(abund2, group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", year == "2023", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

nonoridtTT23$tukey <-      c("a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a")

nonoridtTT23$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

nonoriplotTT23 <-ggplot(nonoridtTT23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Toftrees 2023")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
nonoriplotTT23

nonoridataTT23 <- filter(abund2,  group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", year == "2023", method == "Soil Cores")

zerononoriTT23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=nonoridataTT23, family = nbinom2)
summary(zerononoriTT23)
lsmeans(zerononoriTT23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerononoriTT23, plot = T)
testZeroInflation(zerononoriTT23)

## NON-ORIBATIDA TT 2024

nonoridtTT24 <-filter(abund2, group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", year == "2024", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

nonoridtTT24$tukey <- c("a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a")

nonoridtTT24$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

nonoriplotTT24 <-ggplot(nonoridtTT24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Toftrees 2024")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
nonoriplotTT24

nonoridataTT24 <- filter(abund2, group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", year == "2024", method == "Soil Cores")

zerononoriTT24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=nonoridataTT24, family = nbinom2)
summary(zerononoriTT24)
lsmeans(zerononoriTT24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerononoriTT24, plot = T)
testZeroInflation(zerononoriTT24)

## NON-ORIBATIDA MV 2023

nonoridtMV23 <-filter(abund2, group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", year == "2023", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

nonoridtMV23$tukey <- c("a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","ab","a","a","a",
                       "a","a","a","a","a",
                       "a","b","a","a","a",
                       "a","a","a","a","a")

nonoridtMV23$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

nonoriplotMV23 <-ggplot(nonoridtMV23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Mountain View 2023")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
nonoriplotMV23

nonoridataMV23 <- filter(abund2,group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", year == "2023", method == "Soil Cores")

zerononoriMV23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=nonoridataMV23, family = nbinom2)
summary(zerononoriMV23)
lsmeans(zerononoriMV23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerononoriMV23, plot = T)
testZeroInflation(zerononoriMV23)

## NON-ORIBATIDA MV 2024

nonoridtMV24 <-filter(abund2, group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", year == "2024", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

nonoridtMV24$tukey <- c("a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a",
                       "a","a","a","a","a")
nonoridtMV24$tukeytest  <- c("1","2","3","4","5",
                            "6","7","8","9","10",
                            "11","12","13","14","15",
                            "16","17","18","19","20",
                            "21","22","23","24","25",
                            "26","27","28","29","30")

nonoriplotMV24 <-ggplot(nonoridtMV24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Mountain View 2024")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
nonoriplotMV24

nonoridataMV24 <- filter(abund2, group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", year == "2024", method == "Soil Cores")

zerononoriMV24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=nonoridataMV24, family = nbinom2)
summary(zerononoriMV24)
lsmeans(zerononoriMV24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerononoriMV24, plot = T)
testZeroInflation(zerononoriMV24)

###NON-ORIBATIDA TT 2023+2024

nonori2yrdtTT <- filter(abund2, group == "Non-oribatid Mites" | group == "Predatory Mites", site == "TT", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

nonori2yrdtTT$tukeytest <- c("a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a")
nonori2yrplotTT <- 
  ggplot(nonori2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Non-Oribatid Mite Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
nonori2yrplotTT

###NON ORIBATID MV 2023+2024

nonori2yrdtMV <- filter(abund2, group == "Non-oribatid Mites" | group == "Predatory Mites", site == "MV", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

nonori2yrdtMV$tukeytest <- c("a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "b","a","a","a","a",
                            "a","a","ab","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","b","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a",
                            "a","a","a","a","a")
nonori2yrdtMV$tukey <- c("1","2","3","4","5",
                        "6","7","8","9","10",
  "11","12","13","14","15",
  "16","17","18","19","20",
  "21","22","23","24","25",
  "26","27","28","29","30",
  "31","32","33","34","35",
  "36","37","38","39","40",
  "41","42","43","44","45",
  "46","47","48","49","50",
  "51","52","53","54","55",
  "56","57","58","59","60")
nonori2yrplotMV <- 
  ggplot(nonori2yrdtMV, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Non-Oribatid Mites Abundance: Mountain View 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
nonori2yrplotMV

### ORIBATID MITES TT 2023
oridtTT23 <-filter(abund2, group == "Oribatid Mites", site == "TT", year == "2023", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

oridtTT23$tukey <-      c("a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a")

oridtTT23$tukeytest  <- c("1","2","3","4","5",
                             "6","7","8","9","10",
                             "11","12","13","14","15",
                             "16","17","18","19","20",
                             "21","22","23","24","25",
                             "26","27","28","29","30")

oriplotTT23 <-ggplot(oridtTT23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Toftrees 2023")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
oriplotTT23

oridataTT23 <- filter(abund2,  group == "Oribatid Mites", site == "TT", year == "2023", method == "Soil Cores")

zerooriTT23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=oridataTT23, family = nbinom2)
summary(zerooriTT23)
lsmeans(zerooriTT23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerooriTT23, plot = T)
testZeroInflation(zerooriTT23)

## ORIBATIDA TT 2024

oridtTT24 <-filter(abund2, group == "Oribatid Mites", site == "TT", year == "2024", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

oridtTT24$tukey <- c("a","a","b","a","a",
                        "a","a","ab","a","a",
                        "a","a","ab","a","a",
                        "a","a","ab","a","a",
                        "a","a","b","a","a",
                        "a","a","a","a","a")

oridtTT24$tukeytest  <- c("1","2","3","4","5",
                             "6","7","8","9","10",
                             "11","12","13","14","15",
                             "16","17","18","19","20",
                             "21","22","23","24","25",
                             "26","27","28","29","30")

oriplotTT24 <-ggplot(oridtTT24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Toftrees 2024")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
oriplotTT24

oridataTT24 <- filter(abund2, group == "Oribatid Mites", site == "TT", year == "2024", method == "Soil Cores")

zerooriTT24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=oridataTT24, family = nbinom2)
summary(zerooriTT24)
lsmeans(zerooriTT24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerooriTT24, plot = T)
testZeroInflation(zerooriTT24)

## ORIBATIDA MV 2023

oridtMV23 <-filter(abund2, group == "Oribatid Mites", site == "MV", year == "2023", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

oridtMV23$tukey <- c("a","a","a","a","a",
                        "a","a","a","a","a",
                        "a","a","a","a","a",
                        "a","a","a","a","a",
                        "a","a","a","a","a",
                        "a","a","a","a","a")

oridtMV23$tukeytest  <- c("1","2","3","4","5",
                             "6","7","8","9","10",
                             "11","12","13","14","15",
                             "16","17","18","19","20",
                             "21","22","23","24","25",
                             "26","27","28","29","30")

oriplotMV23 <-ggplot(oridtMV23, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Mountain View 2023")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
oriplotMV23

oridataMV23 <- filter(abund2,group == "Oribatid Mites", site == "MV", year == "2023", method == "Soil Cores")

zerooriMV23 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=oridataMV23, family = nbinom2)
summary(zerooriMV23)
lsmeans(zerooriMV23, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerooriMV23, plot = T)
testZeroInflation(zerooriMV23)

## ORIBATIDA MV 2024

oridtMV24 <-filter(abund2, group == "Oribatid Mites", site == "MV", year == "2024", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year, sum) %>%
  group_by(treatment, time.rep) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

oridtMV24$tukey <- c("a","b","a","a","a",
                        "a","b","a","a","a",
                        "a","ab","a","a","a",
                        "a","a","a","a","a",
                        "a","b","a","a","a",
                        "a","ab","a","a","a")
oridtMV24$tukeytest  <- c("1","2","3","4","5",
                             "6","7","8","9","10",
                             "11","12","13","14","15",
                             "16","17","18","19","20",
                             "21","22","23","24","25",
                             "26","27","28","29","30")

oriplotMV24 <-ggplot(oridtMV24, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Acari Abundance: Mountain View 2024")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukey), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")
oriplotMV24

oridataMV24 <- filter(abund2, group == "Oribatid Mites", site == "MV", year == "2024", method == "Soil Cores")

zerooriMV24 <- glmmTMB(count ~ treatment + time.rep + treatment*time.rep + (1|time.rep), data=oridataMV24, family = nbinom2)
summary(zerooriMV24)
lsmeans(zerooriMV24, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerooriMV24, plot = T)
testZeroInflation(zerooriMV24)

###ORIBATIDA TT 2023+2024

ori2yrdtTT <- filter(abund2, group == "Oribatid Mites", site == "TT", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

ori2yrdtTT$tukeytest <-    c("a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "ab","a","a","a","a",
                             "a","a","a","a","a",
                             "ab","a","a","a","a",
                             "a","a","a","a","a",
                             "ab","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "b","a","a","a","a")

ori2yrdtTT$tukey <-    c("1","2","3","4","5",
                         "6","7","8","9","10",
                         "11","12","13","14","15",
                         "16","17","18","19","20",
                         "21","22","23","24","25",
                         "26","27","28","29","30",
                         "31","32","33","34","35",
                         "36","37","38","39","40",
                         "41","42","43","44","45",
                         "46","47","48","49","50",
                         "51","52","53","54","55",
                         "56","57","58","59","60")
ori2yrplotTT <- 
  ggplot(ori2yrdtTT, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Oribatid Mite Abundance: Toftrees 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
ori2yrplotTT

###ORIBATID MITES MV 2023+2024

ori2yrdtMV <- filter(abund2, group == "Oribatid Mites", site == "MV", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + year + site, sum) %>%
  group_by(treatment, time.rep, year, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

ori2yrdtMV$tukeytest <-   c("a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","ab","a",
                             "a","a","a","a","a",
                             "a","a","a","b","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","a","a",
                             "a","a","a","ab","a",
                             "a","a","a","a","a")

ori2yrdtMV$tukey <-    c("1","2","3","4","5",
                         "6","7","8","9","10",
                         "11","12","13","14","15",
                         "16","17","18","19","20",
                         "21","22","23","24","25",
                         "26","27","28","29","30",
                         "31","32","33","34","35",
                         "36","37","38","39","40",
                         "41","42","43","44","45",
                         "46","47","48","49","50",
                         "51","52","53","54","55",
                         "56","57","58","59","60")
ori2yrplotMV <- 
  ggplot(ori2yrdtMV, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Oribatid Mites Abundance: Mountain View 2023+2024")+
  theme_minimal()+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  facet_grid( rows = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  theme(panel.grid = element_line(colour = "gray65"))
ori2yrplotMV

































###NON-ORIBATID MITES BOTH YEARS COMBINED SITES SEPARATED

nonorbdata2yr <- filter(abund2, group == "Non-Oribatid Mites" , method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

nonorbdata2yr$tukeytest <-  c("a","a","a","a","a",
                              "a","a","a","b","a",
                              "a","a","a","a","a",
                              "a","a","a","ab","a",
                              "a","a","a","a","a",
                              "a","a","a","ab","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","ab","a",
                              "a","a","a","a","a",
                              "a","a","a","ab","a")
nonorbplot2yr <- 
  ggplot(nonorbdata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Non-Oribatid Acari Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.15, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
nonorbplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
nonorbdata2yrTT <- filter(abund2, group == "Non-Oribatid Mites", method == "Soil Cores", site == "TT")

zerononorb2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=nonorbdata2yrTT, family = nbinom2)
summary(zerononorb2yrTT)
lsmeans(zerononorb2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerononorb2yrTT, plot = T)
testZeroInflation(zerononorb2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
nonorbdata2yrMV <- filter(abund2, group == "Non-Oribatid Mites" , method == "Soil Cores", site == "MV")

zerononorb2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=nonorbdata2yrMV, family = nbinom2)
summary(zerononorb2yrMV)
lsmeans(zerononorb2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerononorb2yrMV, plot = T)
testZeroInflation(zerononorb2yrMV)

###ORIBATID MITES BOTH YEARS COMBINED SITES SEPARATED

orbdata2yr <- filter(abund2, group == "Oribatid Mites", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

orbdata2yr$tukeytest <-  c("a","a","3","a","a",
                              "a","a","a","ab","ab",
                              "a","a","a","a","a",
                              "a","a","a","a","b",
                              "a","a","a","a","a",
                              "a","a","a","ab","ab",
                              "a","a","a","a","a",
                              "a","a","a","ab","ab",
                              "a","a","a","a","a",
                              "a","a","a","ab","ab",
                              "a","a","a","a","a",
                              "a","a","a","b","a")
orbplot2yr <- 
  ggplot(orbdata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Oribatid Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.55, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
orbplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
orbdata2yrTT <- filter(abund2, group == "Oribatid Mites", method == "Soil Cores", site == "TT")

zeroorb2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=orbdata2yrTT, family = nbinom2)
summary(zeroorb2yrTT)
lsmeans(zeroorb2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroorb2yrTT, plot = T)
testZeroInflation(zeroorb2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
orbdata2yrMV <- filter(abund2, group == "Oribatid Mites", method == "Soil Cores", site == "MV")

zeroorb2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=orbdata2yrMV, family = nbinom2)
summary(zeroorb2yrMV)
lsmeans(zeroorb2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeroorb2yrMV, plot = T)
testZeroInflation(zeroorb2yrMV)


###PREDATORY MITES BOTH YEARS COMBINED SITES SEPARATED

predmitedata2yr <- filter(abund2, group == "Predatory Mites", method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

predmitedata2yr$tukeytest <-  c("a","a","a","a","a",
                                "a","a","a","a","a",
                                "a","a","a","a","a",
                                "a","a","a","a","a",
                                "a","a","a","a","a",
                                "a","a","a","a","a",
                                "a","a","a","a","a",
                                "a","a","a","a","a",
                                "a","a","b","a","a",
                                "a","a","a","a","a",
                                "a","a","a","a","a",
                                "a","a","a","a","a")
predmiteplot2yr <- 
  ggplot(predmitedata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Predatory Acari Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.55, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
predmiteplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
predmitedata2yrTT <- filter(abund2, group == "Predatory Mites", method == "Soil Cores", site == "TT")

zeropredmite2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=predmitedata2yrTT, family = nbinom2)
summary(zeropredmite2yrTT)
lsmeans(zeropredmite2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeropredmite2yrTT, plot = T)
testZeroInflation(zeropredmite2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
predmitedata2yrMV <- filter(abund2, group == "Predatory Mites" , method == "Soil Cores", site == "MV")

zeropredmite2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=predmitedata2yrMV, family = nbinom2)
summary(zeropredmite2yrMV)
lsmeans(zeropredmite2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zeropredmite2yrMV, plot = T)
testZeroInflation(zeropredmite2yrMV)

###SPIDERS BOTH YEARS COMBINED SITES SEPARATED

spiderdata2yr <- filter(abund2, group == "Araneae" , method == "Pitfalls") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

spiderdata2yr$tukeytest <-  c("a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a",
                              "a","a","a","a","a")
spiderplot2yr <- 
  ggplot(spiderdata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Araneae Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.55, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
spiderplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
spiderdata2yrTT <- filter(abund2, group == "Araneae", method == "Pitfalls", site == "TT")

zerospider2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=spiderdata2yrTT, family = nbinom2)
summary(zerospider2yrTT)
lsmeans(zerospider2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerospider2yrTT, plot = T)
testZeroInflation(zerospider2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
spiderdata2yrMV <- filter(abund2, group == "Araneae" , method == "Pitfalls", site == "MV")

zerospider2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=spiderdata2yrMV, family = nbinom2)
summary(zerospider2yrMV)
lsmeans(zerospider2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerospider2yrMV, plot = T)
testZeroInflation(zerospider2yrMV)

###INSECT LARVAE BOTH YEARS COMBINED SITES SEPARATED

larvaedata2yr <- filter(abund2, group == "Insect Larvae" , method == "Pitfalls" | method == "Soil Cores") %>%
  aggregate(count ~ treatment + time.rep + rep.num + site, sum) %>%
  group_by(treatment, time.rep, site) %>%
  dplyr::summarise(w=mean(count), sd=sd(count))

larvaedata2yr$tukeytest <-  c("ab","a","b","a","b",
                              "b","a","b","a","a",
                              "ab","a","ab","a","a",
                              "ab","a","ab","a","a",
                              "a","a","ab","a","a",
                              "a","a","ab","a","a",
                              "ab","a","ab","a","a",
                              "a","a","ab","a","a",
                              "ab","a","b","a","a",
                              "a","a","b","a","a",
                              "b","a","ab","a","a",
                              "a","a","ab","a","a")
larvaeplot2yr <- 
  ggplot(larvaedata2yr, aes(x=time.rep, y=w, fill=treatment))+
  geom_bar(stat="identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Application", y = "Mean Abundance ± SE")+
  ggtitle("Effect of Insecticides on Insect Larvae Abundance: 2023+2024")+
  theme_minimal()+
  theme(legend.position = c(.55, 0.95),
        legend.justification = c("right","top"),
        legend.box.just = "right",
        legend.margin = margin(6,6,6,6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  guides(fill = guide_legend(title = "Treatment"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3,
            vjust = -0.7, hjust = -0.15, colour = "gray25")+
  facet_grid(cols = vars(site))
larvaeplot2yr

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **TOFTREES**
larvaedata2yrTT <- filter(abund2, group == "Insect Larvae", method == "Pitfalls" | method == "Soil Cores", site == "TT")

zerolarvae2yrTT <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=larvaedata2yrTT, family = nbinom2)
summary(zerolarvae2yrTT)
lsmeans(zerolarvae2yrTT, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerolarvae2yrTT, plot = T)
testZeroInflation(zerolarvae2yrTT)

##ZERO INFLATED MODEL USING YEAR AND TIME REP AS RANDOM EFFECTS **MOUNTAIN VIEW**
larvaedata2yrMV <- filter(abund2, group == "Insect Larvae" , method == "Pitfalls" | method == "Soil Cores", site == "MV")

zerolarvae2yrMV <- glmmTMB(count ~ treatment + time.rep + time.rep*treatment + (1|time.rep+year), data=larvaedata2yrMV, family = nbinom2)
summary(zerolarvae2yrMV)
lsmeans(zerolarvae2yrMV, pairwise ~ treatment | time.rep, adjust = "none")

simulateResiduals(fittedModel = zerolarvae2yrMV, plot = T)
testZeroInflation(zerolarvae2yrMV)

###PIE CHART SHOWING ARTHROPOD BREAKDOWN BY SAMPLING###

piechartpitfalldt <- abund2 %>% filter(method == "Pitfalls") %>% filter( group != "Other" , group != "Diplopoda", group != "Insect Larvae", group != "") %>% aggregate(count ~ group + method, sum)

piechartpitfalldt2 <- piechartpitfalldt %>% dplyr::mutate(group = dplyr::recode(group, "Predatory Mites" = "Acari", "Oribatid Mites" = "Acari", "Non-Oribatid Mites" = "Acari", "Carabidae" = "Coleoptera", "Staphylinidae" = "Coleoptera", "Elateridae" = "Coleoptera")) %>% aggregate(count ~ group + method, sum)


piechartpitfalldt3 <- piechartpitfalldt2 %>%
  mutate(perc = round(count/ sum(count), digits = 4)) %>%
  mutate(labels = scales::percent(perc)) %>%
  mutate(group = fct_rev(fct_inorder(group)))%>%
  mutate(text_y = cumsum(count) - count/2)

pitfalls <- ggplot(piechartpitfalldt3, aes(x = "" , y = count, fill = group)) +
  geom_col(width = 1, color = 1, alpha = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ffa600","#ff6e54","#dd5182","#955196","#003f5c"))+
  geom_label_repel(aes(label = labels, y = text_y), 
                   nudge_x = 0.6, nudge_y = 0.6,
                   size = 5, show.legend = F, alpha = 1,max.overlaps = getOption("ggrepel.max.overlaps", default = 10), colour = "white", segment.color = "black")+
  guides(fill = guide_legend(title = "Taxa")) +
  theme_void()+
  ggtitle("Pitfalls (n=168,838)")

pitfalls

 ###SOIL CORES###

piechartsoilcoredt <- abund2 %>% filter(method == "Soil Cores") %>% filter( group != "Other" , group != "Diplopoda", group != "Insect Larvae", group != "") %>% aggregate(count ~ group + method, sum)

piechartsoilcoredt2 <- piechartsoilcoredt %>% dplyr::mutate(group = dplyr::recode(group, "Predatory Mites" = "Acari", "Oribatid Mites" = "Acari", "Non-Oribatid Mites" = "Acari","Carabidae" = "Coleoptera", "Staphylinidae" = "Coleoptera", "Elateridae" = "Coleoptera")) %>% aggregate(count ~ group + method, sum)


piechartsoilcoredt3 <- piechartsoilcoredt2 %>%
  mutate(perc = round(count/ sum(count), digits = 4)) %>%
  mutate(labels = scales::percent(perc)) %>%
  mutate(group = fct_rev(fct_inorder(group)))%>%
  mutate(text_y = cumsum(count) - count/2) 

soilcores <- ggplot(piechartsoilcoredt3, aes(x = "" , y = count, fill = group)) +
  geom_col(width = 1, color = 1, alpha = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#ffa600","#ff6e54","#dd5182","#955196","#003f5c"))+
  geom_label_repel(aes(label = labels, y = text_y), 
                   nudge_x = 0.6, nudge_y = -0.6,
                   size = 5, show.legend = F, alpha = 1, max.overlaps = getOption("ggrepel.max.overlaps", default = 10), colour = "white", segment.color = "black")+
  guides(fill = guide_legend(title = "Taxa")) +
  theme_void()+
  ggtitle("Soil Cores (n=60,217)")
soilcores

soilcoren <- abund2 %>% aggregate(count ~ method, sum)

combinedpie <- grid.arrange(pitfalls, soilcores, ncol = 2)
combinedpie 

####DIVERSITY COMPONENT

div1 <- rowid_to_column(abund2, "ID")

divtrt2 <- div1 %>%
  unite("taxa", lowtaxon:morph, remove = FALSE) %>%
  aggregate(count ~ treatment + taxa + rep.num, sum)

divtrt3 <- divtrt2 %>%
  pivot_wider(names_from = taxa,
              values_from = count,
              values_fn = list) %>%
  mutate_all(~replace(., lengths(.) == 0 , 0))


divtrt3[,4] <- as.numeric(unlist(divtrt3[,4]))
divtrt3[,5] <- as.numeric(unlist(divtrt3[,5]))
divtrt3[,6] <- as.numeric(unlist(divtrt3[,6]))
divtrt3[,7] <- as.numeric(unlist(divtrt3[,7]))
divtrt3[,8] <- as.numeric(unlist(divtrt3[,8]))
divtrt3[,9] <- as.numeric(unlist(divtrt3[,9]))
divtrt3[,10] <- as.numeric(unlist(divtrt3[,10]))
divtrt3[,11] <- as.numeric(unlist(divtrt3[,11]))
divtrt3[,12] <- as.numeric(unlist(divtrt3[,12]))
divtrt3[,13] <- as.numeric(unlist(divtrt3[,13]))
divtrt3[,14] <- as.numeric(unlist(divtrt3[,14]))
divtrt3[,15] <- as.numeric(unlist(divtrt3[,15]))
divtrt3[,16] <- as.numeric(unlist(divtrt3[,16]))
divtrt3[,17] <- as.numeric(unlist(divtrt3[,17]))
divtrt3[,18] <- as.numeric(unlist(divtrt3[,18]))
divtrt3[,19] <- as.numeric(unlist(divtrt3[,19]))
divtrt3[,20] <- as.numeric(unlist(divtrt3[,20]))
divtrt3[,21] <- as.numeric(unlist(divtrt3[,21]))
divtrt3[,22] <- as.numeric(unlist(divtrt3[,22]))
divtrt3[,23] <- as.numeric(unlist(divtrt3[,23]))
divtrt3[,24] <- as.numeric(unlist(divtrt3[,24]))
divtrt3[,25] <- as.numeric(unlist(divtrt3[,25]))
divtrt3[,26] <- as.numeric(unlist(divtrt3[,26]))
divtrt3[,27] <- as.numeric(unlist(divtrt3[,27]))
divtrt3[,28] <- as.numeric(unlist(divtrt3[,28]))
divtrt3[,29] <- as.numeric(unlist(divtrt3[,29]))
divtrt3[,30] <- as.numeric(unlist(divtrt3[,30]))
divtrt3[,31] <- as.numeric(unlist(divtrt3[,31]))
divtrt3[,32] <- as.numeric(unlist(divtrt3[,32]))
divtrt3[,33] <- as.numeric(unlist(divtrt3[,33]))
divtrt3[,34] <- as.numeric(unlist(divtrt3[,34]))
divtrt3[,35] <- as.numeric(unlist(divtrt3[,35]))
divtrt3[,36] <- as.numeric(unlist(divtrt3[,36]))
divtrt3[,37] <- as.numeric(unlist(divtrt3[,37]))
divtrt3[,38] <- as.numeric(unlist(divtrt3[,38]))
divtrt3[,39] <- as.numeric(unlist(divtrt3[,39]))
divtrt3[,40] <- as.numeric(unlist(divtrt3[,40]))
divtrt3[,41] <- as.numeric(unlist(divtrt3[,41]))
divtrt3[,42] <- as.numeric(unlist(divtrt3[,42]))
divtrt3[,43] <- as.numeric(unlist(divtrt3[,43]))
divtrt3[,44] <- as.numeric(unlist(divtrt3[,44]))
divtrt3[,45] <- as.numeric(unlist(divtrt3[,45]))
divtrt3[,46] <- as.numeric(unlist(divtrt3[,46]))
divtrt3[,47] <- as.numeric(unlist(divtrt3[,47]))
divtrt3[,48] <- as.numeric(unlist(divtrt3[,48]))
divtrt3[,49] <- as.numeric(unlist(divtrt3[,49]))
divtrt3[,50] <- as.numeric(unlist(divtrt3[,50]))
divtrt3[,51] <- as.numeric(unlist(divtrt3[,51]))
divtrt3[,52] <- as.numeric(unlist(divtrt3[,52]))
divtrt3[,53] <- as.numeric(unlist(divtrt3[,53]))
divtrt3[,54] <- as.numeric(unlist(divtrt3[,54]))
divtrt3[,55] <- as.numeric(unlist(divtrt3[,55]))
divtrt3[,56] <- as.numeric(unlist(divtrt3[,56]))
divtrt3[,57] <- as.numeric(unlist(divtrt3[,57]))
divtrt3[,58] <- as.numeric(unlist(divtrt3[,58]))
divtrt3[,59] <- as.numeric(unlist(divtrt3[,59]))
divtrt3[,60] <- as.numeric(unlist(divtrt3[,60]))
divtrt3[,61] <- as.numeric(unlist(divtrt3[,61]))
divtrt3[,62] <- as.numeric(unlist(divtrt3[,62]))
divtrt3[,63] <- as.numeric(unlist(divtrt3[,63]))
divtrt3[,64] <- as.numeric(unlist(divtrt3[,64]))
divtrt3[,65] <- as.numeric(unlist(divtrt3[,65]))
divtrt3[,66] <- as.numeric(unlist(divtrt3[,66]))
divtrt3[,67] <- as.numeric(unlist(divtrt3[,67]))
divtrt3[,68] <- as.numeric(unlist(divtrt3[,68]))
divtrt3[,69] <- as.numeric(unlist(divtrt3[,69]))
divtrt3[,70] <- as.numeric(unlist(divtrt3[,70]))
divtrt3[,71] <- as.numeric(unlist(divtrt3[,71]))
divtrt3[,72] <- as.numeric(unlist(divtrt3[,72]))
divtrt3[,73] <- as.numeric(unlist(divtrt3[,73]))
divtrt3[,74] <- as.numeric(unlist(divtrt3[,74]))
divtrt3[,75] <- as.numeric(unlist(divtrt3[,75]))
divtrt3[,76] <- as.numeric(unlist(divtrt3[,76]))
divtrt3[,77] <- as.numeric(unlist(divtrt3[,77]))
divtrt3[,78] <- as.numeric(unlist(divtrt3[,78]))
divtrt3[,79] <- as.numeric(unlist(divtrt3[,79]))
divtrt3[,80] <- as.numeric(unlist(divtrt3[,80]))
divtrt3[,81] <- as.numeric(unlist(divtrt3[,81]))
divtrt3[,82] <- as.numeric(unlist(divtrt3[,82]))
divtrt3[,83] <- as.numeric(unlist(divtrt3[,83]))
divtrt3[,84] <- as.numeric(unlist(divtrt3[,84]))
divtrt3[,85] <- as.numeric(unlist(divtrt3[,85]))
divtrt3[,86] <- as.numeric(unlist(divtrt3[,86]))
divtrt3[,87] <- as.numeric(unlist(divtrt3[,87]))
divtrt3[,88] <- as.numeric(unlist(divtrt3[,88]))
divtrt3[,89] <- as.numeric(unlist(divtrt3[,89]))
divtrt3[,90] <- as.numeric(unlist(divtrt3[,90]))
divtrt3[,91] <- as.numeric(unlist(divtrt3[,91]))
divtrt3[,92] <- as.numeric(unlist(divtrt3[,92]))
divtrt3[,93] <- as.numeric(unlist(divtrt3[,93]))
divtrt3[,94] <- as.numeric(unlist(divtrt3[,94]))
divtrt3[,95] <- as.numeric(unlist(divtrt3[,95]))
divtrt3[,96] <- as.numeric(unlist(divtrt3[,96]))
divtrt3[,97] <- as.numeric(unlist(divtrt3[,97]))
divtrt3[,98] <- as.numeric(unlist(divtrt3[,98]))
divtrt3[,99] <- as.numeric(unlist(divtrt3[,99]))
divtrt3[,100] <- as.numeric(unlist(divtrt3[,100]))
divtrt3[,101] <- as.numeric(unlist(divtrt3[,101]))
divtrt3[,102] <- as.numeric(unlist(divtrt3[,102]))
divtrt3[,103] <- as.numeric(unlist(divtrt3[,103]))
divtrt3[,104] <- as.numeric(unlist(divtrt3[,104]))
divtrt3[,105] <- as.numeric(unlist(divtrt3[,105]))
divtrt3[,106] <- as.numeric(unlist(divtrt3[,106]))
divtrt3[,107] <- as.numeric(unlist(divtrt3[,107]))
divtrt3[,108] <- as.numeric(unlist(divtrt3[,108]))
divtrt3[,109] <- as.numeric(unlist(divtrt3[,109]))
divtrt3[,110] <- as.numeric(unlist(divtrt3[,110]))
divtrt3[,111] <- as.numeric(unlist(divtrt3[,111]))
divtrt3[,112] <- as.numeric(unlist(divtrt3[,112]))
divtrt3[,113] <- as.numeric(unlist(divtrt3[,113]))
divtrt3[,114] <- as.numeric(unlist(divtrt3[,114]))
divtrt3[,115] <- as.numeric(unlist(divtrt3[,115]))
divtrt3[,116] <- as.numeric(unlist(divtrt3[,116]))
divtrt3[,117] <- as.numeric(unlist(divtrt3[,117]))
divtrt3[,118] <- as.numeric(unlist(divtrt3[,118]))
divtrt3[,119] <- as.numeric(unlist(divtrt3[,119]))
divtrt3[,120] <- as.numeric(unlist(divtrt3[,120]))
divtrt3[,121] <- as.numeric(unlist(divtrt3[,121]))
divtrt3[,122] <- as.numeric(unlist(divtrt3[,122]))


sapply(divtrt3[,4:122], class)


permdata <- divtrt3[,4:122]

permMDS <- metaMDS(permdata, distance = "bray", k =2 , autotransform = TRUE, noshare = TRUE)
permMDS
ordiplot(permMDS, type="n")
orditorp(permMDS, display = "species", col = "red", air = 0.01)
orditorp(permMDS, display = "sites", cex = 1.25, air = 0.01)


trtMDS <- metaMDS(divtrt3[,4:122],noshare = TRUE, autotransform = TRUE, k=2, distance = "bray")
trtMDS
plot(trtMDS)


NMDS1 <- trtMDS$points[,1]
NMDS2 <- trtMDS$points[,2]
plot <- cbind(divtrt3, NMDS1, NMDS2)

p <- ggplot(plot, aes(NMDS1, NMDS2, colour = divtrt3$treatment))+
  geom_point(position=position_jitter(.1), aes(shape = divtrt3$treatment))+
  scale_color_manual(name = "Treatment" , values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000","#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000","#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000","#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000","#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000","#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  stat_ellipse(aes(fill=treatment), alpha=.1,type='t',linewidth =1, geom="polygon", show.legend = FALSE)+
  theme_minimal()+
  labs(title="Insecticide Treatment's Effect on Arthropod Community")+
  scale_shape_discrete(name = "Treatment")
p

stressplot(trtMDS)

div.matrix <- as.matrix(divtrt3[,4:122])
div.dist <- vegdist(div.matrix, method = "bray")
div.dist
set.seed(36)
trt.div <- adonis2(div.dist~as.factor(divtrt3$treatment), data = divtrt3, permutations = 9999)
trt.div
disp.div <- betadisper(div.dist, group = divtrt3$treatment)
disp.div
permutest(disp.div)
anova(disp.div)
plot(disp.div, hull = FALSE, ellipse = TRUE)

###Community composition by site###

divsite1 <- rowid_to_column(abund2, "ID")

divsite2 <- divsite1 %>%
  unite("taxa", lowtaxon:morph, remove = FALSE) %>%
  unite("site_year", site:year, remove = FALSE) %>%
  aggregate(count ~ site + taxa + site_year + rep.num, sum)

divsite3 <- divsite2 %>%
  pivot_wider(names_from = taxa,
              values_from = count,
              values_fn = list) %>%
  mutate_all(~replace(., lengths(.) == 0 , 0))

str(divsite3)


divsite3[,4] <- as.numeric(unlist(divsite3[,4]))
divsite3[,5] <- as.numeric(unlist(divsite3[,5]))
divsite3[,6] <- as.numeric(unlist(divsite3[,6]))
divsite3[,7] <- as.numeric(unlist(divsite3[,7]))
divsite3[,8] <- as.numeric(unlist(divsite3[,8]))
divsite3[,9] <- as.numeric(unlist(divsite3[,9]))
divsite3[,10] <- as.numeric(unlist(divsite3[,10]))
divsite3[,11] <- as.numeric(unlist(divsite3[,11]))
divsite3[,12] <- as.numeric(unlist(divsite3[,12]))
divsite3[,13] <- as.numeric(unlist(divsite3[,13]))
divsite3[,14] <- as.numeric(unlist(divsite3[,14]))
divsite3[,15] <- as.numeric(unlist(divsite3[,15]))
divsite3[,16] <- as.numeric(unlist(divsite3[,16]))
divsite3[,17] <- as.numeric(unlist(divsite3[,17]))
divsite3[,18] <- as.numeric(unlist(divsite3[,18]))
divsite3[,19] <- as.numeric(unlist(divsite3[,19]))
divsite3[,20] <- as.numeric(unlist(divsite3[,20]))
divsite3[,21] <- as.numeric(unlist(divsite3[,21]))
divsite3[,22] <- as.numeric(unlist(divsite3[,22]))
divsite3[,23] <- as.numeric(unlist(divsite3[,23]))
divsite3[,24] <- as.numeric(unlist(divsite3[,24]))
divsite3[,25] <- as.numeric(unlist(divsite3[,25]))
divsite3[,26] <- as.numeric(unlist(divsite3[,26]))
divsite3[,27] <- as.numeric(unlist(divsite3[,27]))
divsite3[,28] <- as.numeric(unlist(divsite3[,28]))
divsite3[,29] <- as.numeric(unlist(divsite3[,29]))
divsite3[,30] <- as.numeric(unlist(divsite3[,30]))
divsite3[,31] <- as.numeric(unlist(divsite3[,31]))
divsite3[,32] <- as.numeric(unlist(divsite3[,32]))
divsite3[,33] <- as.numeric(unlist(divsite3[,33]))
divsite3[,34] <- as.numeric(unlist(divsite3[,34]))
divsite3[,35] <- as.numeric(unlist(divsite3[,35]))
divsite3[,36] <- as.numeric(unlist(divsite3[,36]))
divsite3[,37] <- as.numeric(unlist(divsite3[,37]))
divsite3[,38] <- as.numeric(unlist(divsite3[,38]))
divsite3[,39] <- as.numeric(unlist(divsite3[,39]))
divsite3[,40] <- as.numeric(unlist(divsite3[,40]))
divsite3[,41] <- as.numeric(unlist(divsite3[,41]))
divsite3[,42] <- as.numeric(unlist(divsite3[,42]))
divsite3[,43] <- as.numeric(unlist(divsite3[,43]))
divsite3[,44] <- as.numeric(unlist(divsite3[,44]))
divsite3[,45] <- as.numeric(unlist(divsite3[,45]))
divsite3[,46] <- as.numeric(unlist(divsite3[,46]))
divsite3[,47] <- as.numeric(unlist(divsite3[,47]))
divsite3[,48] <- as.numeric(unlist(divsite3[,48]))
divsite3[,49] <- as.numeric(unlist(divsite3[,49]))
divsite3[,50] <- as.numeric(unlist(divsite3[,50]))
divsite3[,51] <- as.numeric(unlist(divsite3[,51]))
divsite3[,52] <- as.numeric(unlist(divsite3[,52]))
divsite3[,53] <- as.numeric(unlist(divsite3[,53]))
divsite3[,54] <- as.numeric(unlist(divsite3[,54]))
divsite3[,55] <- as.numeric(unlist(divsite3[,55]))
divsite3[,56] <- as.numeric(unlist(divsite3[,56]))
divsite3[,57] <- as.numeric(unlist(divsite3[,57]))
divsite3[,58] <- as.numeric(unlist(divsite3[,58]))
divsite3[,59] <- as.numeric(unlist(divsite3[,59]))
divsite3[,60] <- as.numeric(unlist(divsite3[,60]))
divsite3[,61] <- as.numeric(unlist(divsite3[,61]))
divsite3[,62] <- as.numeric(unlist(divsite3[,62]))
divsite3[,63] <- as.numeric(unlist(divsite3[,63]))
divsite3[,64] <- as.numeric(unlist(divsite3[,64]))
divsite3[,65] <- as.numeric(unlist(divsite3[,65]))
divsite3[,66] <- as.numeric(unlist(divsite3[,66]))
divsite3[,67] <- as.numeric(unlist(divsite3[,67]))
divsite3[,68] <- as.numeric(unlist(divsite3[,68]))
divsite3[,69] <- as.numeric(unlist(divsite3[,69]))
divsite3[,70] <- as.numeric(unlist(divsite3[,70]))
divsite3[,71] <- as.numeric(unlist(divsite3[,71]))
divsite3[,72] <- as.numeric(unlist(divsite3[,72]))
divsite3[,73] <- as.numeric(unlist(divsite3[,73]))
divsite3[,74] <- as.numeric(unlist(divsite3[,74]))
divsite3[,75] <- as.numeric(unlist(divsite3[,75]))
divsite3[,76] <- as.numeric(unlist(divsite3[,76]))
divsite3[,77] <- as.numeric(unlist(divsite3[,77]))
divsite3[,78] <- as.numeric(unlist(divsite3[,78]))
divsite3[,79] <- as.numeric(unlist(divsite3[,79]))
divsite3[,80] <- as.numeric(unlist(divsite3[,80]))
divsite3[,81] <- as.numeric(unlist(divsite3[,81]))
divsite3[,82] <- as.numeric(unlist(divsite3[,82]))
divsite3[,83] <- as.numeric(unlist(divsite3[,83]))
divsite3[,84] <- as.numeric(unlist(divsite3[,84]))
divsite3[,85] <- as.numeric(unlist(divsite3[,85]))
divsite3[,86] <- as.numeric(unlist(divsite3[,86]))
divsite3[,87] <- as.numeric(unlist(divsite3[,87]))
divsite3[,88] <- as.numeric(unlist(divsite3[,88]))
divsite3[,89] <- as.numeric(unlist(divsite3[,89]))
divsite3[,90] <- as.numeric(unlist(divsite3[,90]))
divsite3[,91] <- as.numeric(unlist(divsite3[,91]))
divsite3[,92] <- as.numeric(unlist(divsite3[,92]))
divsite3[,93] <- as.numeric(unlist(divsite3[,93]))
divsite3[,94] <- as.numeric(unlist(divsite3[,94]))
divsite3[,95] <- as.numeric(unlist(divsite3[,95]))
divsite3[,96] <- as.numeric(unlist(divsite3[,96]))
divsite3[,97] <- as.numeric(unlist(divsite3[,97]))
divsite3[,98] <- as.numeric(unlist(divsite3[,98]))
divsite3[,99] <- as.numeric(unlist(divsite3[,99]))
divsite3[,100] <- as.numeric(unlist(divsite3[,100]))
divsite3[,101] <- as.numeric(unlist(divsite3[,101]))
divsite3[,102] <- as.numeric(unlist(divsite3[,102]))
divsite3[,103] <- as.numeric(unlist(divsite3[,103]))
divsite3[,104] <- as.numeric(unlist(divsite3[,104]))
divsite3[,105] <- as.numeric(unlist(divsite3[,105]))
divsite3[,106] <- as.numeric(unlist(divsite3[,106]))
divsite3[,107] <- as.numeric(unlist(divsite3[,107]))
divsite3[,108] <- as.numeric(unlist(divsite3[,108]))
divsite3[,109] <- as.numeric(unlist(divsite3[,109]))
divsite3[,110] <- as.numeric(unlist(divsite3[,110]))
divsite3[,111] <- as.numeric(unlist(divsite3[,111]))
divsite3[,112] <- as.numeric(unlist(divsite3[,112]))
divsite3[,113] <- as.numeric(unlist(divsite3[,113]))
divsite3[,114] <- as.numeric(unlist(divsite3[,114]))
divsite3[,115] <- as.numeric(unlist(divsite3[,115]))
divsite3[,116] <- as.numeric(unlist(divsite3[,116]))
divsite3[,117] <- as.numeric(unlist(divsite3[,117]))
divsite3[,118] <- as.numeric(unlist(divsite3[,118]))
divsite3[,119] <- as.numeric(unlist(divsite3[,119]))
divsite3[,120] <- as.numeric(unlist(divsite3[,120]))
divsite3[,121] <- as.numeric(unlist(divsite3[,121]))
divsite3[,122] <- as.numeric(unlist(divsite3[,122]))
divsite3[,123] <- as.numeric(unlist(divsite3[,123]))


sapply(divsite3[,3:123], class)

siteMDS <- metaMDS(divsite3[,4:123],noshare = TRUE, autotransform = TRUE, k=2, distance = "bray")
siteMDS
plot(siteMDS)


NMDS1A <- siteMDS$points[,1]
NMDS2A <- siteMDS$points[,2]
plot2 <- cbind(divsite3, NMDS1A, NMDS2A)


p2 <- ggplot(plot2, aes(NMDS1A, NMDS2A, colour = divsite3$site_year))+
  geom_point(position=position_jitter(.1), aes(shape = divsite3$site_year))+
  stat_ellipse(aes(fill=site), alpha=.15,type='t',linewidth =1, geom="polygon", show.legend = FALSE)+
  theme_minimal()+
  labs(title="Differences in Arthropod Community across Site and Year")+
  scale_shape_discrete(name = "Site/Year Combinations")+
  scale_color_manual(name = "Site/Year Combinations", values = c("#845ec2","#d65db1","#ff9671","#ffc75f"))
p2

stressplot(siteMDS)

div.site.matrix <- as.matrix(divsite3[,4:123])
div.site.dist <- vegdist(div.site.matrix, method = "bray")
div.site.dist
set.seed(36)
site.div <- adonis2(div.site.dist~as.factor(divsite3$site_year), data = divsite3, permutations = 9999)
site.div
disp.site.div <- betadisper(div.site.dist, group = divsite3$site)
disp.site.div
permutest(disp.site.div)
anova(disp.site.div)
plot(disp.site.div, hull = TRUE, ellipse = FALSE)

###community composition by sampling date###

divtime1 <- rowid_to_column(abund2, "ID")

divtime2 <- divtime1 %>%
  unite("taxa", lowtaxon:morph, remove = FALSE) %>%
  aggregate(count ~ time.rep + taxa + rep.num, sum)

divtime3 <- divtime2 %>%
  pivot_wider(names_from = taxa,
              values_from = count,
              values_fn = list) %>%
  mutate_all(~replace(., lengths(.) == 0 , 0))


divtime3[,2] <- as.numeric(unlist(divtime3[,2]))
divtime3[,3] <- as.numeric(unlist(divtime3[,3]))
divtime3[,4] <- as.numeric(unlist(divtime3[,4]))
divtime3[,5] <- as.numeric(unlist(divtime3[,5]))
divtime3[,6] <- as.numeric(unlist(divtime3[,6]))
divtime3[,7] <- as.numeric(unlist(divtime3[,7]))
divtime3[,8] <- as.numeric(unlist(divtime3[,8]))
divtime3[,9] <- as.numeric(unlist(divtime3[,9]))
divtime3[,10] <- as.numeric(unlist(divtime3[,10]))
divtime3[,11] <- as.numeric(unlist(divtime3[,11]))
divtime3[,12] <- as.numeric(unlist(divtime3[,12]))
divtime3[,13] <- as.numeric(unlist(divtime3[,13]))
divtime3[,14] <- as.numeric(unlist(divtime3[,14]))
divtime3[,15] <- as.numeric(unlist(divtime3[,15]))
divtime3[,16] <- as.numeric(unlist(divtime3[,16]))
divtime3[,17] <- as.numeric(unlist(divtime3[,17]))
divtime3[,18] <- as.numeric(unlist(divtime3[,18]))
divtime3[,19] <- as.numeric(unlist(divtime3[,19]))
divtime3[,20] <- as.numeric(unlist(divtime3[,20]))
divtime3[,21] <- as.numeric(unlist(divtime3[,21]))
divtime3[,22] <- as.numeric(unlist(divtime3[,22]))
divtime3[,23] <- as.numeric(unlist(divtime3[,23]))
divtime3[,24] <- as.numeric(unlist(divtime3[,24]))
divtime3[,25] <- as.numeric(unlist(divtime3[,25]))
divtime3[,26] <- as.numeric(unlist(divtime3[,26]))
divtime3[,27] <- as.numeric(unlist(divtime3[,27]))
divtime3[,28] <- as.numeric(unlist(divtime3[,28]))
divtime3[,29] <- as.numeric(unlist(divtime3[,29]))
divtime3[,30] <- as.numeric(unlist(divtime3[,30]))
divtime3[,31] <- as.numeric(unlist(divtime3[,31]))
divtime3[,32] <- as.numeric(unlist(divtime3[,32]))
divtime3[,33] <- as.numeric(unlist(divtime3[,33]))
divtime3[,34] <- as.numeric(unlist(divtime3[,34]))
divtime3[,35] <- as.numeric(unlist(divtime3[,35]))
divtime3[,36] <- as.numeric(unlist(divtime3[,36]))
divtime3[,37] <- as.numeric(unlist(divtime3[,37]))
divtime3[,38] <- as.numeric(unlist(divtime3[,38]))
divtime3[,39] <- as.numeric(unlist(divtime3[,39]))
divtime3[,40] <- as.numeric(unlist(divtime3[,40]))
divtime3[,41] <- as.numeric(unlist(divtime3[,41]))
divtime3[,42] <- as.numeric(unlist(divtime3[,42]))
divtime3[,43] <- as.numeric(unlist(divtime3[,43]))
divtime3[,44] <- as.numeric(unlist(divtime3[,44]))
divtime3[,45] <- as.numeric(unlist(divtime3[,45]))
divtime3[,46] <- as.numeric(unlist(divtime3[,46]))
divtime3[,47] <- as.numeric(unlist(divtime3[,47]))
divtime3[,48] <- as.numeric(unlist(divtime3[,48]))
divtime3[,49] <- as.numeric(unlist(divtime3[,49]))
divtime3[,50] <- as.numeric(unlist(divtime3[,50]))
divtime3[,51] <- as.numeric(unlist(divtime3[,51]))
divtime3[,52] <- as.numeric(unlist(divtime3[,52]))
divtime3[,53] <- as.numeric(unlist(divtime3[,53]))
divtime3[,54] <- as.numeric(unlist(divtime3[,54]))
divtime3[,55] <- as.numeric(unlist(divtime3[,55]))
divtime3[,56] <- as.numeric(unlist(divtime3[,56]))
divtime3[,57] <- as.numeric(unlist(divtime3[,57]))
divtime3[,58] <- as.numeric(unlist(divtime3[,58]))
divtime3[,59] <- as.numeric(unlist(divtime3[,59]))
divtime3[,60] <- as.numeric(unlist(divtime3[,60]))
divtime3[,61] <- as.numeric(unlist(divtime3[,61]))
divtime3[,62] <- as.numeric(unlist(divtime3[,62]))
divtime3[,63] <- as.numeric(unlist(divtime3[,63]))
divtime3[,64] <- as.numeric(unlist(divtime3[,64]))
divtime3[,65] <- as.numeric(unlist(divtime3[,65]))
divtime3[,66] <- as.numeric(unlist(divtime3[,66]))
divtime3[,67] <- as.numeric(unlist(divtime3[,67]))
divtime3[,68] <- as.numeric(unlist(divtime3[,68]))
divtime3[,69] <- as.numeric(unlist(divtime3[,69]))
divtime3[,70] <- as.numeric(unlist(divtime3[,70]))
divtime3[,71] <- as.numeric(unlist(divtime3[,71]))
divtime3[,72] <- as.numeric(unlist(divtime3[,72]))
divtime3[,73] <- as.numeric(unlist(divtime3[,73]))
divtime3[,74] <- as.numeric(unlist(divtime3[,74]))
divtime3[,75] <- as.numeric(unlist(divtime3[,75]))
divtime3[,76] <- as.numeric(unlist(divtime3[,76]))
divtime3[,77] <- as.numeric(unlist(divtime3[,77]))
divtime3[,78] <- as.numeric(unlist(divtime3[,78]))
divtime3[,79] <- as.numeric(unlist(divtime3[,79]))
divtime3[,80] <- as.numeric(unlist(divtime3[,80]))
divtime3[,81] <- as.numeric(unlist(divtime3[,81]))
divtime3[,82] <- as.numeric(unlist(divtime3[,82]))
divtime3[,83] <- as.numeric(unlist(divtime3[,83]))
divtime3[,84] <- as.numeric(unlist(divtime3[,84]))
divtime3[,85] <- as.numeric(unlist(divtime3[,85]))
divtime3[,86] <- as.numeric(unlist(divtime3[,86]))
divtime3[,87] <- as.numeric(unlist(divtime3[,87]))
divtime3[,88] <- as.numeric(unlist(divtime3[,88]))
divtime3[,89] <- as.numeric(unlist(divtime3[,89]))
divtime3[,90] <- as.numeric(unlist(divtime3[,90]))
divtime3[,91] <- as.numeric(unlist(divtime3[,91]))
divtime3[,92] <- as.numeric(unlist(divtime3[,92]))
divtime3[,93] <- as.numeric(unlist(divtime3[,93]))
divtime3[,94] <- as.numeric(unlist(divtime3[,94]))
divtime3[,95] <- as.numeric(unlist(divtime3[,95]))
divtime3[,96] <- as.numeric(unlist(divtime3[,96]))
divtime3[,97] <- as.numeric(unlist(divtime3[,97]))
divtime3[,98] <- as.numeric(unlist(divtime3[,98]))
divtime3[,99] <- as.numeric(unlist(divtime3[,99]))
divtime3[,100] <- as.numeric(unlist(divtime3[,100]))
divtime3[,101] <- as.numeric(unlist(divtime3[,101]))
divtime3[,102] <- as.numeric(unlist(divtime3[,102]))
divtime3[,103] <- as.numeric(unlist(divtime3[,103]))
divtime3[,104] <- as.numeric(unlist(divtime3[,104]))
divtime3[,105] <- as.numeric(unlist(divtime3[,105]))
divtime3[,106] <- as.numeric(unlist(divtime3[,106]))
divtime3[,107] <- as.numeric(unlist(divtime3[,107]))
divtime3[,108] <- as.numeric(unlist(divtime3[,108]))
divtime3[,109] <- as.numeric(unlist(divtime3[,109]))
divtime3[,110] <- as.numeric(unlist(divtime3[,110]))
divtime3[,111] <- as.numeric(unlist(divtime3[,111]))
divtime3[,112] <- as.numeric(unlist(divtime3[,112]))
divtime3[,113] <- as.numeric(unlist(divtime3[,113]))
divtime3[,114] <- as.numeric(unlist(divtime3[,114]))
divtime3[,115] <- as.numeric(unlist(divtime3[,115]))
divtime3[,116] <- as.numeric(unlist(divtime3[,116]))
divtime3[,117] <- as.numeric(unlist(divtime3[,117]))
divtime3[,118] <- as.numeric(unlist(divtime3[,118]))
divtime3[,119] <- as.numeric(unlist(divtime3[,119]))
divtime3[,120] <- as.numeric(unlist(divtime3[,120]))
divtime3[,121] <- as.numeric(unlist(divtime3[,121]))
divtime3[,122] <- as.numeric(unlist(divtime3[,122]))

sapply(divtime3[,2:122], class)

timeMDS <- metaMDS(divtime3[,3:122],noshare = TRUE, autotransform = TRUE, k=2, distance = "bray")
timeMDS
plot(siteMDS)

NMDS1B <- timeMDS$points[,1]
NMDS2B <- timeMDS$points[,2]
plot3 <- cbind(divtime3, NMDS1B, NMDS2B)

p3 <- ggplot(plot3, aes(NMDS1B, NMDS2B, colour = as.factor(divtime3$time.rep)))+
  geom_point(position=position_jitter(.1), shape =2, aes(shape = divetime3$time.rep))+
  stat_ellipse(aes(fill=time.rep), alpha=.15,type='t',linewidth =1, geom="polygon", show.legend = FALSE)+
  theme_minimal()+
  scale_color_manual(name = "Sampling Period" ,values = c("#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  labs(title="Differences in Arthropod Community across Sampling Period")
p3

stressplot(timeMDS)

div.time.matrix <- as.matrix(divtime3[,3:122])
div.time.dist <- vegdist(div.time.matrix, method = "bray")
div.time.dist
set.seed(36)
time.div <- adonis2(div.time.dist~as.factor(divtime3$time.rep), data = divsite3, permutations = 9999)
time.div
disp.time.div <- betadisper(div.time.dist, group = divtime3$time.rep)
disp.time.div
permutest(disp.time.div)
anova(disp.time.div)
plot(disp.time.div, hull = FALSE, ellipse = TRUE)

###SHANNON AND SIMPSON DIVERSITY

divdata1 <- rowid_to_column(abund2, "ID")

divdata2 <- divdata1 %>%
  unite("taxa", lowtaxon:morph, remove = FALSE) %>%
  unite("site_year", site:year, remove = FALSE) %>%
  aggregate(count ~ taxa + site_year + treatment, sum)

divdata3 <- divdata2 %>%
  pivot_wider(names_from = taxa,
              values_from = count,
              values_fn = list) %>%
  mutate_all(~replace(., lengths(.) == 0 , 0))


divdata3[,3] <- as.numeric(unlist(divdata3[,3]))
divdata3[,4] <- as.numeric(unlist(divdata3[,4]))
divdata3[,5] <- as.numeric(unlist(divdata3[,5]))
divdata3[,6] <- as.numeric(unlist(divdata3[,6]))
divdata3[,7] <- as.numeric(unlist(divdata3[,7]))
divdata3[,8] <- as.numeric(unlist(divdata3[,8]))
divdata3[,9] <- as.numeric(unlist(divdata3[,9]))
divdata3[,10] <- as.numeric(unlist(divdata3[,10]))
divdata3[,11] <- as.numeric(unlist(divdata3[,11]))
divdata3[,12] <- as.numeric(unlist(divdata3[,12]))
divdata3[,13] <- as.numeric(unlist(divdata3[,13]))
divdata3[,14] <- as.numeric(unlist(divdata3[,14]))
divdata3[,15] <- as.numeric(unlist(divdata3[,15]))
divdata3[,16] <- as.numeric(unlist(divdata3[,16]))
divdata3[,17] <- as.numeric(unlist(divdata3[,17]))
divdata3[,18] <- as.numeric(unlist(divdata3[,18]))
divdata3[,19] <- as.numeric(unlist(divdata3[,19]))
divdata3[,20] <- as.numeric(unlist(divdata3[,20]))
divdata3[,21] <- as.numeric(unlist(divdata3[,21]))
divdata3[,22] <- as.numeric(unlist(divdata3[,22]))
divdata3[,23] <- as.numeric(unlist(divdata3[,23]))
divdata3[,24] <- as.numeric(unlist(divdata3[,24]))
divdata3[,25] <- as.numeric(unlist(divdata3[,25]))
divdata3[,26] <- as.numeric(unlist(divdata3[,26]))
divdata3[,27] <- as.numeric(unlist(divdata3[,27]))
divdata3[,28] <- as.numeric(unlist(divdata3[,28]))
divdata3[,29] <- as.numeric(unlist(divdata3[,29]))
divdata3[,30] <- as.numeric(unlist(divdata3[,30]))
divdata3[,31] <- as.numeric(unlist(divdata3[,31]))
divdata3[,32] <- as.numeric(unlist(divdata3[,32]))
divdata3[,33] <- as.numeric(unlist(divdata3[,33]))
divdata3[,34] <- as.numeric(unlist(divdata3[,34]))
divdata3[,35] <- as.numeric(unlist(divdata3[,35]))
divdata3[,36] <- as.numeric(unlist(divdata3[,36]))
divdata3[,37] <- as.numeric(unlist(divdata3[,37]))
divdata3[,38] <- as.numeric(unlist(divdata3[,38]))
divdata3[,39] <- as.numeric(unlist(divdata3[,39]))
divdata3[,40] <- as.numeric(unlist(divdata3[,40]))
divdata3[,41] <- as.numeric(unlist(divdata3[,41]))
divdata3[,42] <- as.numeric(unlist(divdata3[,42]))
divdata3[,43] <- as.numeric(unlist(divdata3[,43]))
divdata3[,44] <- as.numeric(unlist(divdata3[,44]))
divdata3[,45] <- as.numeric(unlist(divdata3[,45]))
divdata3[,46] <- as.numeric(unlist(divdata3[,46]))
divdata3[,47] <- as.numeric(unlist(divdata3[,47]))
divdata3[,48] <- as.numeric(unlist(divdata3[,48]))
divdata3[,49] <- as.numeric(unlist(divdata3[,49]))
divdata3[,50] <- as.numeric(unlist(divdata3[,50]))
divdata3[,51] <- as.numeric(unlist(divdata3[,51]))
divdata3[,52] <- as.numeric(unlist(divdata3[,52]))
divdata3[,53] <- as.numeric(unlist(divdata3[,53]))
divdata3[,54] <- as.numeric(unlist(divdata3[,54]))
divdata3[,55] <- as.numeric(unlist(divdata3[,55]))
divdata3[,56] <- as.numeric(unlist(divdata3[,56]))
divdata3[,57] <- as.numeric(unlist(divdata3[,57]))
divdata3[,58] <- as.numeric(unlist(divdata3[,58]))
divdata3[,59] <- as.numeric(unlist(divdata3[,59]))
divdata3[,60] <- as.numeric(unlist(divdata3[,60]))
divdata3[,61] <- as.numeric(unlist(divdata3[,61]))
divdata3[,62] <- as.numeric(unlist(divdata3[,62]))
divdata3[,63] <- as.numeric(unlist(divdata3[,63]))
divdata3[,64] <- as.numeric(unlist(divdata3[,64]))
divdata3[,65] <- as.numeric(unlist(divdata3[,65]))
divdata3[,66] <- as.numeric(unlist(divdata3[,66]))
divdata3[,67] <- as.numeric(unlist(divdata3[,67]))
divdata3[,68] <- as.numeric(unlist(divdata3[,68]))
divdata3[,69] <- as.numeric(unlist(divdata3[,69]))
divdata3[,70] <- as.numeric(unlist(divdata3[,70]))
divdata3[,71] <- as.numeric(unlist(divdata3[,71]))
divdata3[,72] <- as.numeric(unlist(divdata3[,72]))
divdata3[,73] <- as.numeric(unlist(divdata3[,73]))
divdata3[,74] <- as.numeric(unlist(divdata3[,74]))
divdata3[,75] <- as.numeric(unlist(divdata3[,75]))
divdata3[,76] <- as.numeric(unlist(divdata3[,76]))
divdata3[,77] <- as.numeric(unlist(divdata3[,77]))
divdata3[,78] <- as.numeric(unlist(divdata3[,78]))
divdata3[,79] <- as.numeric(unlist(divdata3[,79]))
divdata3[,80] <- as.numeric(unlist(divdata3[,80]))
divdata3[,81] <- as.numeric(unlist(divdata3[,81]))
divdata3[,82] <- as.numeric(unlist(divdata3[,82]))
divdata3[,83] <- as.numeric(unlist(divdata3[,83]))
divdata3[,84] <- as.numeric(unlist(divdata3[,84]))
divdata3[,85] <- as.numeric(unlist(divdata3[,85]))
divdata3[,86] <- as.numeric(unlist(divdata3[,86]))
divdata3[,87] <- as.numeric(unlist(divdata3[,87]))
divdata3[,88] <- as.numeric(unlist(divdata3[,88]))
divdata3[,89] <- as.numeric(unlist(divdata3[,89]))
divdata3[,90] <- as.numeric(unlist(divdata3[,90]))
divdata3[,91] <- as.numeric(unlist(divdata3[,91]))
divdata3[,92] <- as.numeric(unlist(divdata3[,92]))
divdata3[,93] <- as.numeric(unlist(divdata3[,93]))
divdata3[,94] <- as.numeric(unlist(divdata3[,94]))
divdata3[,95] <- as.numeric(unlist(divdata3[,95]))
divdata3[,96] <- as.numeric(unlist(divdata3[,96]))
divdata3[,97] <- as.numeric(unlist(divdata3[,97]))
divdata3[,98] <- as.numeric(unlist(divdata3[,98]))
divdata3[,99] <- as.numeric(unlist(divdata3[,99]))
divdata3[,100] <- as.numeric(unlist(divdata3[,100]))
divdata3[,101] <- as.numeric(unlist(divdata3[,101]))
divdata3[,102] <- as.numeric(unlist(divdata3[,102]))
divdata3[,103] <- as.numeric(unlist(divdata3[,103]))
divdata3[,104] <- as.numeric(unlist(divdata3[,104]))
divdata3[,105] <- as.numeric(unlist(divdata3[,105]))
divdata3[,106] <- as.numeric(unlist(divdata3[,106]))
divdata3[,107] <- as.numeric(unlist(divdata3[,107]))
divdata3[,108] <- as.numeric(unlist(divdata3[,108]))
divdata3[,109] <- as.numeric(unlist(divdata3[,109]))
divdata3[,110] <- as.numeric(unlist(divdata3[,110]))
divdata3[,111] <- as.numeric(unlist(divdata3[,111]))
divdata3[,112] <- as.numeric(unlist(divdata3[,112]))
divdata3[,113] <- as.numeric(unlist(divdata3[,113]))
divdata3[,114] <- as.numeric(unlist(divdata3[,114]))
divdata3[,115] <- as.numeric(unlist(divdata3[,115]))
divdata3[,116] <- as.numeric(unlist(divdata3[,116]))
divdata3[,117] <- as.numeric(unlist(divdata3[,117]))
divdata3[,118] <- as.numeric(unlist(divdata3[,118]))
divdata3[,119] <- as.numeric(unlist(divdata3[,119]))
divdata3[,120] <- as.numeric(unlist(divdata3[,120]))
divdata3[,121] <- as.numeric(unlist(divdata3[,121]))
divdata3[,122] <- as.numeric(unlist(divdata3[,122]))



divdataMV24 <- filter(divdata3, site_year == "MV_2024")
  
specnumberMV24 <- specnumber(divdataMV24[,3:122])
specnumberMV24


divdataMV23 <- filter(divdata3, site_year == "MV_2023")

specnumberMV23 <- specnumber(divdataMV23[,3:122])
specnumberMV23


divdataTT23 <- filter(divdata3, site_year == "TT_2023")

specnumberTT23 <- specnumber(divdataTT23[,3:122])
specnumberTT23


divdataTT24 <- filter(divdata3, site_year == "TT_2024")

specnumberTT24 <- specnumber(divdataTT24[,3:122])
specnumberTT24

###SHANNON
shannondivMV24 <- diversity(divdataMV24[,3:122], index = "shannon")
shannondivMV24

shannondivMV23 <- diversity(divdataMV23[,3:122], index = "shannon")
shannondivMV23

shannondivTT24 <- diversity(divdataTT24[,3:122], index = "shannon")
shannondivTT24

shannondivTT23 <- diversity(divdataTT23[,3:122], index = "shannon")
shannondivTT23

#SIMPSON

simpsondivMV24 <- diversity(divdataMV24[,3:122], index = "simpson")
simpsondivMV24

simpsondivMV23 <- diversity(divdataMV23[,3:122], index = "simpson")
simpsondivMV23

simpsondivTT24 <- diversity(divdataTT24[,3:122], index = "simpson")
simpsondivTT24

simpsondivTT23 <- diversity(divdataTT23[,3:122], index = "simpson")
simpsondivTT23

###BOOTSTRAPPING CONFIDENCE INTERVALS FOR SIMPSON AND SHANNON DIVERSITY INDEXES FOR EACH TREATMENT FOR EACH SITE_YEAR COMBINATION

divdata4 <- divdata3 %>%
  unite("site_year_trt", site_year:treatment, remove = TRUE) 

bootstrap1 <- sbdiv(divdata4[,2:100], f = divdata4$site_year_trt, theta = c("Shannon"), method = c("asht"), type = "Sequen" )
bootstrap1

bootstrap2 <- sbdiv(divdata4[,2:100], f = divdata4$site_year_trt, theta = c("Simpson"), method = c("asht"), type = "Sequen" )
bootstrap2

#obtaining overall arthropod counts for each site and year
nMV23 <- filter(abund2, year == "2023", site == "MV") %>% aggregate(count ~ treatment, sum)
nMV24 <- filter(abund2, year == "2024", site == "MV") %>% aggregate(count ~ treatment, sum)
nTT23 <- filter(abund2, year == "2023", site == "TT") %>% aggregate(count ~ treatment, sum)
nTT24 <- filter(abund2, year == "2024", site == "TT") %>% aggregate(count ~ treatment, sum)


