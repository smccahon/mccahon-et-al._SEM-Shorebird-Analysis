#----------------------------------------#
#           SEM Model Building           #
#      Spring 2023 Invertebrate Dataset  #
#  Created by Shelby McCahon on 12/10/25 #
#          Modified on 05/25/2026       #
#----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(DHARMa)
library(AICcmodavg)

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

invert <- read.csv("Data/Invert_2023.csv")

# subset to spring (n = 35)
invert <- subset(invert, Season == "Spring")

#------------------------------------------------------------------------------#
#                        convert factors to numeric                         ----                        
#------------------------------------------------------------------------------# 

invert <- invert %>% 
  mutate(WaterNeonicDetection = case_when(
    WaterNeonicDetection == "Y" ~ 1,
    WaterNeonicDetection == "N" ~ 0,
    TRUE ~ NA_real_),
    InvertPesticideDetection = case_when(
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

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 

cor_mat <- cor(invert[sapply(invert, is.numeric)], 
               use = "pairwise.complete.obs")

# Convert matrix to a data frame of pairs
cor_df <- as.data.frame(as.table(cor_mat))

# Rename columns for clarity
names(cor_df) <- c("Var1", "Var2", "Correlation")

# Convert factors to characters
cor_df$Var1 <- as.character(cor_df$Var1)
cor_df$Var2 <- as.character(cor_df$Var2)

# Keep only unique pairs (upper triangle)
cor_df <- cor_df[cor_df$Var1 < cor_df$Var2, ]

# Filter by threshold
high_corr <- subset(cor_df, abs(Correlation) > 0.6)

# Sort by correlation strength
high_corr <- high_corr[order(-abs(high_corr$Correlation)), ]

high_corr

#                         Var1                     Var2 Correlation
# 286                 TDS_mg.L             WaterQuality       1.000
# 288       Conductivity_uS.cm             WaterQuality       1.000
# 287             Salinity_ppt             WaterQuality       1.000
# 528       Conductivity_uS.cm                 TDS_mg.L       1.000
# 527             Salinity_ppt                 TDS_mg.L       0.999
# 552       Conductivity_uS.cm             Salinity_ppt       0.999
# 131                 Buffered    NearestCropDistance_m       0.851
# 117             EnvDetection     WaterNeonicDetection       0.793
# 96        Conductivity_uS.cm     PesticideInvert_ng.g       0.779
# 268     PesticideInvert_ng.g             WaterQuality       0.772
# 532     PesticideInvert_ng.g             Salinity_ppt       0.772
# 508     PesticideInvert_ng.g                 TDS_mg.L       0.766
# 325                Diversity InvertPesticideDetection       0.758
# 222    NearestCropDistance_m                PercentAg      -0.739
# 227                 Buffered                PercentAg      -0.733
# 291                  Biomass                Diversity       0.717
# 350 InvertPesticideDetection               Permanence       0.683
# 110 InvertPesticideDetection     WaterNeonicDetection      -0.667

#------------------------------------------------------------------------------#
#                      which variables can go in model                      ----                        
#------------------------------------------------------------------------------# 

table(invert$Buffered) # 27 not buffered, 8 buffered -> just use % ag
table(invert$WaterNeonicDetection) # 14 with detections, 21 without detections
table(invert$InvertPesticideDetection) # 4 with detections, 6 without

table(invert$EnvDetection) # 18 with detections 
# this includes the 4 wetlands with invert pesticide detections 
# (that didn't have neonics in the wetland)

table(invert$Permanence) # very uneven, not worth including

# variables of interest and for simplicity:
# % cropland, biomass, diversity, water neonic conc

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# add a small constant to allow gamma distribution use
# I chose the smallest constant I could without causing model to fail
# 0.001 is what I used for other invert analysis
# piecewiseSEM does not support tweedie distribution which is ideal dist.
# next best is gamma
# log also works, but you get similar results anyways, and gamma fits data
# better

invert <- invert %>% 
  mutate(Biomass.m = Biomass + 0.001)

# subset to positive values only (n = 14)
# add a small constant so we can log transform
invert.pos <- invert %>% 
  filter(NeonicWater_ng.L > 0) %>% 
  mutate(LogWaterNeonic.m = log(NeonicWater_ng.L + 0.00001))

# specify model
m1 <- glm(Biomass.m ~ LogWaterNeonic.m + PercentAg,
          data = invert.pos,
          family = Gamma(link = "log"))

# # extract standardized coefficients manually from gamma distribution
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.169 is standardized estimate

beta <- coef(m1)["LogWaterNeonic.m"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.147 is standardized estimate

# ...diversity model ----

# use poisson > negative binomial for consistency with other invert analysis
# dispersion = 1.30
m2 <- glm(Diversity ~ LogWaterNeonic.m + PercentAg,
          data = invert.pos,
          family = poisson())

# plot(m2) # seems fine

# ...water neonic conc. model ----

# linear model seems fine
m3 <- lm(LogWaterNeonic.m ~ PercentAg,
         data = invert.pos)

# plot(m3) # seems fine

# run piecewiseSEM
model <- psem(m1,m2,m3)
summary(model, conserve = TRUE)

#------------------------------------------------------------------------------#
#                         VIF of component models                           ----                        
#------------------------------------------------------------------------------# 

# both under 3
vif(m1)
vif(m2)

  