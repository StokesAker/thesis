##ANT MOUNDING

install.packages("multcomp")
library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(multcomp)
library(emmeans)
library(lme4)
library(tidyverse)
library(plyr)
library(lsmeans)

#read in file
ants <- read.csv(file = "https://raw.githubusercontent.com/StokesAker/thesis/main/AntMoundingMaster_v2.csv")

#refactoring variables
ants$rep.num <- as.factor(ants$rep.num)
ants$time.rep <- as.factor(ants$time.rep)

#renaming variables to make it look better
ants1 <- ants %>%
  mutate(treatment = dplyr::recode(treatment, 
                          "control" = "Control",
                          "nov" = "Novaluron",
                          "bif" = "Bifenthrin",
                          "combo" = "Combination")) %>%
  mutate(time.rep = dplyr::recode(time.rep,
                           "0" = "Pre-treatment",
                           "1" = "3 DAT",
                           "2" = "7 DAT",
                           "3" = "10 DAT",
                           "4" = "14 DAT",
                           "5" = "21 DAT")) 

antsdt <- ants1 %>% 
  aggregate(mound.count ~ treatment + time.rep + rep.num + year, sum) %>% 
  group_by(treatment, time.rep, year) %>%
  dplyr::summarise(w=mean(mound.count), sd = sd(mound.count))

antsdt$treatment <- factor(antsdt$treatment, levels = c("Control","Novaluron","Bifenthrin","Combination"))


antsdt$ tukeytest <- c("a","a","b","b","b",
                       "b","bc","b","bc","bc",
                       "a","ab","a","a","b",
                       "b","b","b","c","c",
                       "c","c","a","b","a",
                       "a","a","a","a","a",
                       "a","a","a","ab","a",
                       "a","a","a","a","a",
                       "a","a","ab","a","ab",
                       "a","a","a")



antsplot <-
  ggplot(antsdt, aes(x=time.rep, y= w, fill = treatment))+
  geom_bar(stat= "identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin = w, ymax = w+sd), position = position_dodge(0.9), width = 0.25)+
  labs(x="Days After Treatment", y = "Mean Ant Mounds ± SD", color = "Treatment")+
  ggtitle("Effect of Insecticide Treatment on Ant Mounding")+
  theme_minimal()+
  guides(fill = guide_legend(title = "Treatment"))+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FFA500","#008000"))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 4,
            vjust = -0.8, hjust = -0.5, colour = "gray25")+
  facet_grid(rows = vars(year))+
  theme(panel.grid = element_line(colour = "gray65"))
antsplot


### 2023 ANOVA SQUARE ROOT TRANSFORMATION

antsdata1 <- filter(ants1, year == "2023")

antsdata1$sqrt.mound.count <- sqrt(antsdata1$mound.count)

antanova23 <- aov(sqrt.mound.count ~ treatment*time.rep, data=antsdata1)
summary(antanova23)

plot(antanova23, which = 2)
hist(antanova23$residuals)
shapiro.test(antanova23$residuals)
plot(antanova23, which = 3)


antsemm23 <- emmeans(antanova23, ~treatment*time.rep)
pairs(antsemm23, by = "time.rep", adjust = "none")
pairs(antsemm23, by = "treatment", adjust = "none")

### 2024 ANOVA SQUARE ROOT TRANSFORMATION

antsdata3 <- filter(ants1, year == "2024")

antsdata3$sqrt.mound.count <- sqrt(antsdata3$mound.count)


antanova24 <- aov(sqrt.mound.count ~ treatment*time.rep, data=antsdata3)
summary(antanova24)

plot(antanova24, which = 2)
hist(antanova24$residuals)
shapiro.test(antanova24$residuals)
plot(antanova24, which = 3)


antsemm24 <- emmeans(antanova24, ~treatment*time.rep)
pairs(antsemm24, by = "time.rep", adjust = "none")
pairs(antsemm24, by = "treatment", adjust = "none")

###Sentinel Prey Study


install.packages("caret")
install.packages("pROC")
install.packages("lsmeans")
install.packages("lme4")
library(lsmeans)
library(pROC)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(multcomp)
library(emmeans)
library(lme4)
library(tidyverse)
library(plyr)
library(car)

sentinel <- read.csv(file = "https://raw.githubusercontent.com/StokesAker/thesis/main/SentinelMaster_Final_v6.csv")


sentinel$time.rep <- as.factor(sentinel$time.rep)


sentinel2 <- sentinel %>% 
  dplyr::mutate(treatment = dplyr::recode(treatment,
                                          control = "Control",
                                          nov = "Novaluron",
                                          ace = "Acephate",
                                          bif = "Bifenthrin",
                                          combo = "Combination",
                                          cloth = "Clothianidin")) %>%
  mutate(time.rep = dplyr::recode(time.rep,
                                  "0" = "Pre-treatment",
                                  "1" = "0 DAT",
                                  "2" = "7 DAT",
                                  "3" = "14 DAT",
                                  "4" = "28 DAT",
                                  "5" = "56 DAT")) 

sentinel2$treatment <- as.factor(sentinel2$treatment)

str(sentinel2)

#changing the levels of the treatment factor so everything is compared to controls later on
sentinel2$treatment <- factor(sentinel2$treatment, levels = c("Control","Novaluron","Acephate","Bifenthrin","Clothianidin","Combination"))

levels(sentinel2$treatment)


#I want to merge partial and whole predation into a binary value for predation in general
sentinel2$binary.pred <- ifelse(sentinel2$pred.lvl>0,1,0)


#creating a data frame that gives % predation at the plot level


sentinel3 <- sentinel2 %>% filter(binary.pred == "0" , time.rep != "0 DAT") %>%
  group_by(treatment, time.rep,rep.num,binary.pred, year) %>% dplyr::summarise(count=(n()/5)) 

sentinel3$ppred <- 1-sentinel3$count



#taking the average of predation across replicates for each time point
sentinel4 <- sentinel3 %>%
  dplyr::group_by(treatment, time.rep, year, .drop=FALSE) %>%
  dplyr::summarise(mean = mean(ppred), sd = sd(ppred))


#graphing in ggplot

sentinel5 <- dplyr::filter(sentinel4, time.rep != "0 DAT")

sentinel5$tukeytest <- c("b","a","a","ac","a",
                         "b","ab","ac","a","a",
                         "a","a","a","acd","ac",
                         "ab","ab","ac","a","a",
                         "a","a","a","b","abc",
                         "ab","b","a","a","a",
                         "ab","a","a","abcd","bc",
                         "ab","a","b","a","a",
                         "ab","a","a","bd","abc",
                         "a","ab","bc","a","a",
                         "ab","a","a","c","b",
                         "ab","a","a","a","a")

sentbarplot <- ggplot(sentinel5, aes(fill= treatment, y=mean, x=time.rep))+
  geom_bar(stat = "identity", position = "dodge", alpha = 0.95, colour = "black")+
  geom_errorbar(aes(ymin=mean, ymax = sd+mean), position=position_dodge(0.9), width = 0.25)+
  labs(x="Days After Treatment", y= "Mean Proportion of Predation ± SD")+
  ggtitle("Effect of Insecticide Treatment on Proportion of Sentinel Prey Consumed")+
  theme_minimal()+
  guides(fill=guide_legend(title="Treatment"))+
  scale_fill_manual(values = c("#d3d3d3","#000000","#FF0000","#FFA500","#4169e1","#008000"))+
  facet_grid(row = vars(year))+
  geom_text(aes(label=tukeytest), position = position_dodge(.9), size = 3.5,
            vjust = -0.5, hjust = -.15, colour = "gray1")+
  theme(panel.grid = element_line(colour = "gray65"))

sentbarplot

#ARCSIN TRANSFORM 2023 DATA

sentinel23 <- filter(sentinel3, year == "2023", time.rep != "0 DAT")

sentinel23$asin.ppred <- asin(sqrt(sentinel23$ppred))
sentanova5 <- aov(asin.ppred ~ treatment*time.rep, data=sentinel23)
summary(sentanova5)

plot(sentanova5, which = 2)
hist(sentanova5$residuals)
shapiro.test(sentanova5$residuals)
plot(sentanova5 , which = 3)
leveneTest(sentanova5)

sentemm5 <- emmeans(sentanova5, ~treatment*time.rep)
pairs(sentemm5, by = "time.rep", adjust = "none")
pairs(sentemm5, by = "treatment", adjust = "none")

#ARCSIN TRANSFORM 2024 DATA

sentinel24 <- filter(sentinel3, year == "2024", time.rep != "0 DAT")

sentinel24$asin.ppred <- asin(sqrt(sentinel24$ppred))
sentanova6 <- aov(asin.ppred ~ treatment*time.rep, data=sentinel24)
summary(sentanova6)

plot(sentanova6, which = 2)
hist(sentanova6$residuals)
shapiro.test(sentanova6$residuals)
plot(sentanova6 , which = 3)
leveneTest(sentanova6)

sentemm6 <- emmeans(sentanova6, ~treatment*time.rep)
pairs(sentemm6, by = "time.rep", adjust = "none")
pairs(sentemm6, by = "treatment", adjust = "none")

