#----------------------------------------#
#            SEM Model Building          #
#         Invertebrate 2023 Dataset      #
#  Created by Shelby McCahon on 12/10/25 #
#         Modified on 05/26/2026         #
#----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(DHARMa)
library(AICcmodavg)
library(car)

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

invert <- read.csv("Data/Invert_2023.csv")

#------------------------------------------------------------------------------#
#                        convert factors to numeric                         ----                        
#------------------------------------------------------------------------------# 

invert <- invert %>% 
  mutate(InvertPesticideDetection = case_when(
      InvertPesticideDetection == "Y" ~ 1,
      InvertPesticideDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    EnvDetection = case_when(
      EnvDetection == "Y" ~ 1,
      EnvDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Season = case_when(
      Season == "Spring" ~ 1,
      Season == "Fall" ~ 0),
    Buffered = case_when(
      Buffered == "Y" ~ 1,
      Buffered == "N" ~ 0,
      TRUE ~ NA_real_),
    Permanence = case_when(
      Permanence == "Temporary" ~ 1,
      Permanence == "Seasonal" ~ 2,
      Permanence == "Semipermanent" ~ 3,
      Permanence == "Permanent" ~ 4,
      TRUE ~ NA_real_))

# notes from previous analyses and code:
# I ran a separate exercise to determine which water quality metric to use and 
# I found that pH and conductivity were the most important biologically and
# statistically. I had to log transform conductivity to satisfy model
# assumptions. Correlation code is found in previous scripts but the pasted
# correlations are listed below

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 

# taken from previous code

#                      Var1                     Var2 Correlation
# 552    Conductivity_uS.cm             Salinity_ppt       0.999
# 287          Salinity_ppt             WaterQuality       0.998
# 288    Conductivity_uS.cm             WaterQuality       0.997
# 286              TDS_mg.L             WaterQuality       0.991
# 527          Salinity_ppt                 TDS_mg.L       0.981
# 528    Conductivity_uS.cm                 TDS_mg.L       0.979
# 9                  Julian                   Season      -0.970
# 333          EnvDetection InvertPesticideDetection       0.854
# 227              Buffered                PercentAg      -0.748
# 131              Buffered    NearestCropDistance_m       0.731
# 222 NearestCropDistance_m                PercentAg      -0.725
# 15             Permanence                   Season      -0.630
# 13              Diversity                   Season      -0.626
# 345                Julian               Permanence       0.615

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# notes from previous code and analyses:
# The best distribution would be a tweedie distribution but piecewiseSEM
# does not support that. The next best distribution would be gamma but I 
# had to add a small constant (smallest one I could without causing
# model to fail) to use it. I also had to log transform conductivity to 
# satisfy model assumptions. I wanted to test vegetation cover, but 
# ultimately had to remove to avoid model convergence issues.

invert <- invert %>% 
  mutate(Biomass.m = Biomass + 0.001) # n = 77

invert <- invert %>% 
  mutate(LogConductivity = log(Conductivity_uS.cm))

m1 <- glm(Biomass.m ~ PercentAg + Season + EnvDetection + 
            LogConductivity + Permanence + pH_probe +
            Buffered, 
          data = invert,
          family = Gamma(link = "log"))

# extract standardized coefficients manually from gamma distribution
# % surrounding cropland
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.318 is standardized estimate

# buffer presence
beta <- coef(m1)["Buffered"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.081 is standardized estimate

# season
beta <- coef(m1)["Season"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.124 is standardized estimate

# environmental detection
beta <- coef(m1)["EnvDetection"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # 0.0412 is standardized estimate 

# conductivity
beta <- coef(m1)["LogConductivity"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.037 is standardized estimate

# pH
beta <- coef(m1)["pH_probe"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.0134 is standardized estimate

# permanence
beta <- coef(m1)["Permanence"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # 0.047 is standardized estimate

# ...diversity model ----

# notes from previous analyses:
# data is not dispersed enough to need a negative binomial distributions --
# I get warnings.

m2 <- glm(Diversity ~ PercentAg + Season + LogConductivity + 
            EnvDetection + Permanence + pH_probe + Buffered, 
          data = invert,
          family = poisson())

# ...pesticide detection model ----

# notes:
# I considered adding precipitation and snowfall, but it gave spurious
# results so I removed because I wasn't confident in them

m3 <- glm(EnvDetection ~ PercentAg + Season + Buffered + Permanence,
          data = invert, family = binomial())

# ...water quality index model ----

m4 <- lm(LogConductivity ~ PercentAg + Buffered + Season + Permanence, 
        data = invert)

# ...pH model ----
m5 <- lm(pH_probe ~ PercentAg + Buffered + Season + Permanence, 
         data = invert)

# run piecewiseSEM
model <- psem(m1,m2,m3,m4,m5)
summary(model, conserve = TRUE)

# add correlated errors due to known correlations
model2 <- update(model, 
                 Permanence %~~% Season,
                 PercentAg %~~% Buffered)
summary(model2, conserve = TRUE)

#------------------------------------------------------------------------------#
#                          VIF for component models                         ----                        
#------------------------------------------------------------------------------# 

# all below 3
vif(m1)
vif(m2)
vif(m3)
vif(m4)
vif(m5)

#------------------------------------------------------------------------------#
#                             model diagnostics                             ----                        
#------------------------------------------------------------------------------# 

# m1 ---  good, no severe violations
simulationOutput <- simulateResiduals(fittedModel = m1) # no violations
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # not significant

plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m1)$EnvDetection)  # good
plotResiduals(simulationOutput, form = model.frame(m1)$Season) # significant (maybe due to correlation)
plotResiduals(simulationOutput, form = model.frame(m1)$Permanence) # good
plotResiduals(simulationOutput, form = model.frame(m1)$Conductivity_uS.cm) # good
plotResiduals(simulationOutput, form = model.frame(m1)$pH_probe) # good
plotResiduals(simulationOutput, form = model.frame(m1)$Buffered) # good

# m2 --- all good, no violations
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # not significant

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$Season) # significant (maybe due to correlation)
plotResiduals(simulationOutput, form = model.frame(m2)$Permanence) # good
plotResiduals(simulationOutput, form = model.frame(m2)$Conductivity_uS.cm) # good
plotResiduals(simulationOutput, form = model.frame(m2)$pH_probe) # good
plotResiduals(simulationOutput, form = model.frame(m2)$Buffered) # good

# m3 --- all good, no violations
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m3)$Season) 
plotResiduals(simulationOutput, form = model.frame(m3)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m3)$Permanence) # slight issue

# m4 ---  good, no severe violations with log transformation of response
# patterns, but n.s.
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # pattern but n.s.

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m4)$Season) 
plotResiduals(simulationOutput, form = model.frame(m4)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m4)$Permanence) # pattern but n.s.

# m5 --- good, no severe violations
# patterns in residuals but not significant except for permanence
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m5)$Season) 
plotResiduals(simulationOutput, form = model.frame(m5)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m5)$Permanence) # sign.

