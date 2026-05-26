#-----------------------------------------#
#           SEM Model Building            #
#       Shorebird 2021-2023 Dataset       #
# Created by Shelby McCahon on 12/10/2025 #
#         Modified on 05/26/2026          #
#-----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(DHARMa)
library(AICcmodavg)
library(statmod)
library(MuMIn)
library(semEff)
library(car)

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

birds <- read.csv("Data/Shorebird_2021-2023.csv")

#------------------------------------------------------------------------------#
#                        convert factors to numeric                         ----                        
#------------------------------------------------------------------------------# 

birds <- birds %>% 
  mutate(PlasmaDetection = case_when(
    PlasmaDetection == "Y" ~ 1,
    PlasmaDetection == "N" ~ 0,
    TRUE ~ NA_real_),
    Fat.Binomial = case_when(
      Fat >= 2 ~ 1,
      Fat < 2 ~ 0,
      TRUE ~ NA_real_),
    EnvDetection = case_when(
      EnvDetection == "Y" ~ 1,
      EnvDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Season = case_when(
      Season == "Spring" ~ 1,
      Season == "Fall" ~ 0,
      TRUE ~ NA_real_),
    MigStatus = case_when(
      MigStatus == "Migratory" ~ 1,
      MigStatus == "Resident" ~ 0,
      TRUE ~ NA_real_))

birds <- birds %>%
  filter(complete.cases(PlasmaDetection,
                        Fat, BCI,
                        FatteningIndex,
                        Uric))

# notes from previous code and analyses:
# I originally reduced the dataset down to species with >3 individuals so that
# I could consider species as a random effect. However, the only model that it
# would converge in is the plasma neonic detection model. The fixed effects
# on plasma detection were very similar and in the same direction, but 
# migratory status was significant only in the fixed effects model. However,
# I ultimately decided to remove species as a random effect for consistency 
# with the other bird and invert analysis in 2023 and so that I could include
# the additional 4 birds that I had to remove originally.
# I also ran all of these analyses with Julian and Season, and season performed
# better according to AIC.
# I also tested for correlations in previous scripts and found that the only
# correlated variables of interest were season and drought

# extract one row per site to avoid pseudoreplication in analysis (n = 24 wetlands)
site_data <- birds %>%
  distinct(Site, SPEI, PercentAg, Season, EnvDetection, AnnualSnowfall_in,
           DaysSinceLastPrecipitation_5mm, PrecipitationAmount_7days)

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...wetland pesticide detection model ----

# notes from previous analysis:
# I removed season to avoid overfitting the data. The model without season is 
# also preferred (AICcwt = 0.84).

m1 <- glm(EnvDetection ~ PercentAg + SPEI, 
          data = site_data, 
          family = "binomial")

# ...plasma detection model ----

# notes from previous analyses:
# julian is a better fit than season (AICc wt = 0.81), but I used season
# for consistency with other analyses. Julian has high
# correlation with Season anyways (r = 0.95). Random effect of species 
# worsened model fit (AICc wt in model without random effect = 0.59) so I 
# decided ultimately to remove it. I also considered precipitation variables
# but it did not have a significant or interesting effect.

m2 <- glm(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
                MigStatus + time_hours,
              data = birds, family = "binomial")

# ...body condition index model ----

m3 <- lm(BCI ~ time_hours + PlasmaDetection + SPEI + PercentAg + 
           Season,
         data = birds)


# ...fattening index model ----

m4 <- lm(FatteningIndex ~ MigStatus + Season + SPEI + 
           PercentAg + time_hours + PlasmaDetection + EnvDetection, 
         data = birds)

# ...fat model ----

# notes from previous analysis:
# I had a lot of issues at first with fat as a response variable. I tried
# a lot of different transformations and distributions and ultimately 
# decided to treat fat as binomial (low fat: 0/1 [n = 57] and 
# high fat: 2-5 [n = 27])

m5 <- glm(Fat.Binomial ~ PercentAg + Season + SPEI + PlasmaDetection + 
            MigStatus + time_hours,
          data = birds,
          family = binomial(link = "logit"))

#...uric acid level model ----

m6 <- lm(Uric ~ time_hours + PercentAg + Season + SPEI +
           PlasmaDetection + MigStatus, data = birds)

#------------------------------------------------------------------------------#
#                           run piecewise SEM                               ----                        
#------------------------------------------------------------------------------# 

# run piecewiseSEM
model <- psem(m1,m2,m3,m4,m5,m6)
summary(model, conserve = TRUE)
sem <- update(model,Season %~~% SPEI,
              FatteningIndex %~~% BCI,
              Fat.Binomial %~~% BCI,
              Fat.Binomial %~~% FatteningIndex)

summary(sem, conserve = TRUE)

#------------------------------------------------------------------------------#
#                         VIF of component models                           ----                        
#------------------------------------------------------------------------------# 

# all under 3
vif(m1)
vif(m2)
vif(m3)
vif(m4)
vif(m5)
vif(m6)

#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# wetland pesticide detection model 
# (DHARMa diagnostics unreliable with effect size of 16)
plot(residuals(m1, type = "deviance") ~ fitted(m1))
abline(h = 0, col = "red") # reasonable fit

# plasma detection model --- no severe violations
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m2)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m2)$SPEI) # some pattern

# body condition index model --- all good, no severe patterns
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PlasmaDetection) # significant but no apparent problems
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m3)$SPEI) # good
plotResiduals(simulationOutput, form = model.frame(m3)$FatteningIndex) # good

# fattening index model --- no severe issues
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) # significant
testQuantiles(simulationOutput) # great

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$MigStatus)  # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m4)$BCI) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # good

plot(residuals(m4, type = "deviance") ~ fitted(m4))
abline(h = 0, col = "red") # residuals look just fine though

# fat model --- no severe violations
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # significant

plotResiduals(simulationOutput, form = model.frame(m5)$PercentAg) #  some pattern
plotResiduals(simulationOutput, form = model.frame(m5)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m5)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m5)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m5)$SPEI) # some pattern

# pectoral muscle model --- all good, no patterns
simulationOutput <- simulateResiduals(fittedModel = m6) 
plot(simulationOutput)
testDispersion(m6) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m6)$PercentAg)
plotResiduals(simulationOutput, form = model.frame(m6)$PlasmaDetection)
plotResiduals(simulationOutput, form = model.frame(m6)$time_hours)  
plotResiduals(simulationOutput, form = model.frame(m6)$SPEI)

#------------------------------------------------------------------------------#
#                               odds ratios                                 ----                        
#------------------------------------------------------------------------------# 

# migrant status on plasma detection
exp(coef(m2)) # OR = 0.2646102; migrants have lower odds of detection than residents

1/exp(coef(m2)) # OR = 3.77; residents have higher odds of detection than migrants

# SPEI on plasma detection
exp(coef(m2)) # OR = 2.25; 1 unit increase in SPEI increased odds of detection by 2.25 


# 95% CI for migrant status on detection
1/exp(cbind(Odds_Ratio = coef(m2), confint(m2)))

# 95% CI for SPEI on detection
exp(cbind(Odds_Ratio = coef(m2), confint(m2)))


# SPEI on fat levels
exp(coef(m5)) # OR = 3.04

# 95% CI for SPEI on fat
exp(cbind(Odds_Ratio = coef(m5), confint(m5)))






