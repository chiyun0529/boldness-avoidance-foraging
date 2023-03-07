library(tidyverse)
library(nlme)
library(emmeans)
library(effectsize)
library(merDeriv)
library(outliers)
library(DHARMa)
library(r2glmm)
library(mgcv)
library(GGally)
source("HighstatLibV10.R")

# 1. Boldness: variation and repeatability
boldness<-read_csv("boldness.csv")

# Selecting individuals from the focal populations
boldness2<-subset(boldness,subset = boldness$id>=59)

# Excluding juveniles
boldness_adults<-subset(boldness2,subset = boldness2$sex!="j")

# Testing for multicollinearity among predictors
corvif(boldness_adults[,c(2:4)])

# Variance inflation factors
# 
# GVIF
# sex        1.023728
# population 1.023802
# exp        1.000278

var<-varIdent(form = ~1|as.factor(population))
boldness_glmm<-lme(pc1~as.factor(exp)+as.factor(population)+as.factor(sex),random = ~1|id,data = boldness_adults)
boldness_varIdent<-lme(pc1~as.factor(exp)+as.factor(population)+as.factor(sex),random = ~1|id,weights = varIdent(form = ~1|as.factor(population)),data = boldness_adults)

AIC(boldness_glmm,boldness_varIdent)

# df      AIC
# boldness_glmm      6 357.5860
# boldness_varIdent  7 333.7794

# Model validation
E2<-resid(boldness_varIdent,type = "normalized")
fit<-fitted(boldness_varIdent)

par(mfrow=c(2,3))

plot(x = fit,y = E2)
hist(E2)
plot(as.factor(boldness_adults$exp),E2)
plot(as.factor(boldness_adults$sex),E2)
plot(as.factor(boldness_adults$population),E2)

# Model summary
summary(boldness_varIdent)

# Semi-partial R^2 for the population effect
r2beta(boldness_varIdent,method = "sgv",partial = FALSE)

# Effect   Rsq upper.CL lower.CL
# 1  Model 0.411    0.557    0.279

# 2. Analyses involving avoidance learning
dff5b<-read_csv("dff5b.csv")

# 2.1 First prey attacked
# 2.1.1 GLM using sex, treatment, population, exp, and boldness as predictors
df_1stprey_alltrials<-dff5b[,c(2,3,4,14,15,16)]

prey1_alltrials_start<-glm(`1st_prey2`~.^2,data = df_1stprey_alltrials,family = "binomial")

selection<-step(prey1_alltrials_start,scope = .~.^3,direction = "both")
summary(selection)

prey1_alltrials_best<-glm(`1st_prey2` ~ as.factor(sex) + as.factor(group) + as.factor(exp) + as.factor(population) + 
                            pc1 + as.factor(sex):as.factor(group) + as.factor(sex):as.factor(exp) + as.factor(group):as.factor(exp) + as.factor(exp):as.factor(population) + 
                            as.factor(exp):pc1 + as.factor(population):pc1 + as.factor(exp):as.factor(population):pc1, family = "binomial", 
                          data = df_1stprey_alltrials)

# 2.1.2 Model validation - all good
plot(E<-simulateResiduals(prey1_alltrials_best))

par(mfrow=c(2,3))
plotResiduals(E,as.factor(df_1stprey_alltrials$sex))
plotResiduals(E,as.factor(df_1stprey_alltrials$group))
plotResiduals(E,as.factor(df_1stprey_alltrials$exp))
plotResiduals(E,as.factor(df_1stprey_alltrials$population))
plotResiduals(E,df_1stprey_alltrials$pc1)

# 2.1.3 Summary
summary(prey1_alltrials_best)

# 2.1.4 Effect size
params_glm_prey1<-parameters::model_parameters(prey1_alltrials_best) 
d_prey1<-oddsratio_to_d(params_glm_prey1$Coefficient,log = TRUE)
interpret_prey1<-interpret_cohens_d(d_prey1,rules = "gignac2016")
effectsize_prey1<-data.frame(cbind(params_glm_prey1,d_prey1,interpret_prey1))

write_csv(effectsize_prey1,"effectsize_1stprey_alltrials.csv")

# 2.2 Foraging priority
# 2.2.1 GLMM with ID as a random-effect factor
df_unp50_alltrials<-dff5b[,c(1,2,3,4,8,9,14,15)]
unp50alltrials_glmm<-glmer(cbind(unpalatable_50,total_50-unpalatable_50)~(sex+population+group)*pc1+(1|id),data = df_unp50_alltrials,family = binomial(link = "cloglog"))
summary(unp50alltrials_glmm)
# sigular fit from trivial random effect: variance = 3.61e-14

# 2.2.2 GLM with sex, treatment, population, boldness, and exp as predictors
df_unp50_alltrials<-dff5b[,c(2,3,4,8,9,14,15)]
unp50alltrials_start<-glm(cbind(unpalatable_50,total_50-unpalatable_50)~.^2,data = df_unp50_alltrials,family = binomial(link = "cloglog"))

selection<-step(unp50alltrials_start,scope = .~.^3)
summary(selection)

unp50alltrials_best<-glm(cbind(unpalatable_50, total_50 - unpalatable_50) ~ 
                           as.factor(sex) + as.factor(group) + as.factor(exp) + as.factor(population) + pc1 + as.factor(sex):as.factor(group) + as.factor(sex):as.factor(exp) + 
                           as.factor(sex):as.factor(population) + as.factor(sex):pc1 + as.factor(group):as.factor(exp) + as.factor(group):as.factor(population) + 
                           as.factor(exp):as.factor(population) + as.factor(exp):pc1 + as.factor(population):pc1 + as.factor(exp):as.factor(population):pc1 + 
                           as.factor(sex):as.factor(group):as.factor(population) + as.factor(sex):as.factor(group):as.factor(exp) + as.factor(sex):as.factor(exp):as.factor(population), 
                         family = binomial(link = "cloglog"), data = df_unp50_alltrials)

summary(unp50alltrials_best)

# Model validation - underdispersion/overfitting detected
plot(E<-simulateResiduals(unp50alltrials_best))

par(mfrow=c(2,3))
plotResiduals(E,as.factor(df_unp50_alltrials$sex))
plotResiduals(E,as.factor(df_unp50_alltrials$group))
plotResiduals(E,as.factor(df_unp50_alltrials$exp))
plotResiduals(E,df_unp50_alltrials$pc1)

# 2.2.3 GLM with sex, treatment, boldness, and exp as predictors
df_unp50_alltrials<-dff5b[,c(2,3,4,8,9,15)]
unp50alltrials_start<-glm(cbind(unpalatable_50,total_50-unpalatable_50)~(sex+group+exp)*pc1,data = df_unp50_alltrials,family = binomial(link = "cloglog"))

selection<-step(unp50alltrials_start,scope = .~.^3,direction = "both")
summary(selection)

unp50alltrials_best<-glm(cbind(unpalatable_50, total_50 - unpalatable_50) ~ 
                           as.factor(sex) + as.factor(group) + as.factor(exp) + pc1 + as.factor(sex):pc1, family = binomial(link = "cloglog"), 
                         data = df_unp50_alltrials)

# 2.2.4 Model validation - all good
plot(E<-simulateResiduals(unp50alltrials_best))

par(mfrow=c(2,3))
plotResiduals(E,as.factor(df_unp50_alltrials$sex))
plotResiduals(E,as.factor(df_unp50_alltrials$group))
plotResiduals(E,as.factor(df_unp50_alltrials$exp))
plotResiduals(E,df_unp50_alltrials$pc1)

# 2.2.5 Summary
summary(unp50alltrials_best)

# 2.2.6 Effect size
params_glm_unp50<-parameters::model_parameters(unp50alltrials_best) 
d_unp50<-oddsratio_to_d(params_glm_unp50$Coefficient,log = TRUE)
interpret_unp50<-interpret_cohens_d(d_unp50,rules = "gignac2016")
effectsize_unp50<-data.frame(cbind(params_glm_unp50,d_unp50,interpret_unp50))

write_csv(effectsize_unp50,"effectsize_unp50_alltrials.csv")

# 2.3 Latency
# 2.3.1 Zero-inflated gamma regression with sex, treatment, population, exp, and boldness^2 as predictors
df_zeroif<-dff5b[,c(2,3,4,6,14,15)]

latency_alltrials_start<-glmmTMB(latency~(sex+group+exp+population)*poly(pc1,2),data = df_zeroif,ziformula = ~.,family = ziGamma(link = "log"))

selection<-step(latency_alltrials_start,direction = "both")
summary(selection)
# AIC = 921.2

# Zero-inflated gamma regression with sex, treatment, population, exp, and boldness as predictors
latency_alltrials_start2<-glmmTMB(latency~(sex+group+exp+population)*pc1,data = df_zeroif,ziformula = ~.,family = ziGamma(link = "log"))

selection2<-step(latency_alltrials_start2,direction = "both")
summary(selection2)
# AIC = 918.6, superior

latency_alltrials_best<-glmmTMB(latency ~ as.factor(group) + as.factor(population) + pc1 + as.factor(group):pc1,data = df_zeroif,ziformula = ~.,family = ziGamma(link = "log"))

# 2.3.2 Model validation - all good
plot(E<-simulateResiduals(latency_alltrials_best))

par(mfrow=c(2,3))
plotResiduals(E,as.factor(df_zeroif$group))
plotResiduals(E,as.factor(df_zeroif$population))
plotResiduals(E,df_zeroif$pc1)

# 2.3.3 Summary
summary(latency_alltrials_best)

# 2.4 Overall avoidance
df_unpcons_alltrials<-dff5b[,c(2,3,4,11,13,14,15)]

# 2.4.1 GLM with sex, treatment, population, exp, and boldness as predictors
unpcons_alltrials_start<-glm(cbind(unpalatable_consumed,palatable_consumed)~.^2,data = df_unpcons_alltrials,family = binomial(link = "cloglog"))

selection<-step(unpcons_alltrials_start,scope = .~.^3,direction = "both")

unpcons_alltrials_best<-glm(cbind(unpalatable_consumed, palatable_consumed) ~ 
                              as.factor(sex) + as.factor(group) + as.factor(exp) + as.factor(population) + pc1 + as.factor(sex):as.factor(exp) + as.factor(sex):as.factor(population) + 
                              as.factor(sex):pc1 + as.factor(exp):pc1 + as.factor(sex):as.factor(exp):pc1, family = binomial(link = "cloglog"), 
                            data = df_unpcons_alltrials)


# 2.4.2 Model validation - slight underdispersion b/c unequal variances between sexes
plot(E<-simulateResiduals(unpcons_alltrials_best))

par(mfrow=c(2,3))
plotResiduals(E,as.factor(df_unpcons_alltrials$sex))
plotResiduals(E,as.factor(df_unpcons_alltrials$group))
plotResiduals(E,as.factor(df_unpcons_alltrials$exp))
plotResiduals(E,as.factor(df_unpcons_alltrials$population))
plotResiduals(E,df_unpcons_alltrials$pc1)

# 2.4.3 Summary
summary(unpcons_alltrials_best)

# 2.4.4 Effect size
params_glm_unpcons<-parameters::model_parameters(unpcons_alltrials_best) 
d_unpcons<-oddsratio_to_d(params_glm_unpcons$Coefficient,log = TRUE)
interpret_unpcons<-interpret_cohens_d(d_unpcons,rules = "gignac2016")
effectsize_unpcons<-data.frame(cbind(params_glm_unpcons,d_unpcons,interpret_unpcons))

write_csv(effectsize_unpcons,"effectsize_unpcons_alltrials.csv")