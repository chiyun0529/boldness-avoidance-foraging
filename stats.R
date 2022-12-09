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
dff4b3<-read_csv("dff4b3.csv")

# 2.1 1st_prey
dff4b3.vif<-dff4b3
dff4b3.vif$sex2<-ifelse(dff4b3.vif$sex=="f",1,2)
dff4b3.vif$group2<-ifelse(dff4b3.vif$group=="R",1,2)
dff4b3.vif$pop2<-ifelse(dff4b3.vif$population=="Chai",1,2)

ggpairs(dff4b3.vif[,c(2:4,14,15)],ggplot2::aes(alpha = 0.9, colour = as.factor(group)))
corvif(dff4b3.vif[,c(4,15,19:21)])

# Variance inflation factors
# 
# GVIF
# exp    1.000000
# pc1    1.592233
# sex2   1.040004
# group2 1.015016
# pop2   1.624450

df_1stprey_glmm<-dff4b3[,c(1,2,3,4,14,15,16)]

# This is just to show the triviality of random effect
reg_prey1_glmm<-glmer(`1st_prey2`~.^2+(id|exp),data = df_1stprey_glmm,family = "binomial")
VarCorr(reg_prey1_glmm)
# Random Std. Dev = 3.16e-9

# GLM
df_1stprey<-dff4b3[,c(2,3,4,14,15,16)]

reg_prey1_all<-glm(`1st_prey2`~.^2,data = df_1stprey,family = "binomial")

# Model selection
selection<-step(reg_prey1_all,scope = .~.^3,direction = "both")
summary(selection)

# Best model
reg_prey1_best<-glm(formula = `1st_prey2` ~ sex + group + exp + population + 
                      pc1 + sex:group + sex:population + sex:pc1 + group:exp + 
                      group:pc1 + exp:population + exp:pc1 + population:pc1 + exp:population:pc1 + 
                      sex:population:pc1, family = "binomial", data = df_1stprey)

# Checking for overdispersion
overdisp_fun(reg_prey1_best)

# chisq      ratio        rdf          p 
# 57.2814906  0.9546915 60.0000000  0.5756889 

# Calculating McFadden's pseudo r^2 for logistic regression
with(summary(reg_prey1_best), 1 - deviance/null.deviance)

# Model validation
plot(E<-simulateResiduals(reg_prey1_best))

par(mfrow=c(2,3))
plotResiduals(E,as.factor(df_1stprey$sex))
plotResiduals(E,as.factor(df_1stprey$group))
plotResiduals(E,as.factor(df_1stprey$exp))
plotResiduals(E,as.factor(df_1stprey$population))
plotResiduals(E,df_1stprey$pc1)

# Model summary
summary(reg_prey1_best)

# Effect size
params_reg_prey1<-parameters::model_parameters(reg_prey1_best) 
d_prey1<-oddsratio_to_d(params_reg_prey1$Coefficient,log = TRUE)
interpret_prey1<-interpret_cohens_d(d_prey1,rules = "gignac2016")
effectsize_prey1<-data.frame(cbind(params_reg_prey1,d_prey1,interpret_prey1))


# 2.2 Foraging priority

# Setting up data frame for GLM
df_unp50<-dff4b3[,c(2,3,4,8,10,14,15)]

reg_unp50_all<-glm(cbind(unpalatable_50,palatable_50)~.^2,data = df_unp50,family = binomial(link = "cloglog"))

# Model selection
selection<-step(reg_unp50_all,scope = .~.^3)
summary(selection)

# Best model
reg_unp50_best<-glm(formula = cbind(unpalatable_50, palatable_50) ~ sex + group + 
                      exp + population + pc1 + sex:pc1 + group:pc1 + exp:population + 
                      exp:pc1 + population:pc1 + exp:population:pc1, family = binomial(link = "cloglog"), 
                    data = df_unp50)

# Checking for overdispersion
overdisp_fun(reg_unp50_best)

# chisq      ratio        rdf          p 
# 65.6084992  1.0251328 64.0000000  0.4208006 

# McFadden's pseudo R^2
with(summary(reg_unp50_best), 1 - deviance/null.deviance)

# Model validation
plot(E<-simulateResiduals(reg_unp50_best))

par(mfrow=c(2,3))
plotResiduals(E,as.factor(df_unp50$sex))
plotResiduals(E,as.factor(df_unp50$group))
plotResiduals(E,as.factor(df_unp50$exp))
plotResiduals(E,as.factor(df_unp50$population))
plotResiduals(E,df_unp50$pc1)

# Model summary
summary(reg_unp50_best)

# Effect size
params_reg_unp50<-parameters::model_parameters(reg_unp50_best) 
d_unp50<-oddsratio_to_d(params_reg_unp50$Coefficient,log = TRUE)
interpret_unp50<-interpret_cohens_d(d_unp50,rules = "gignac2016")
effectsize_unp50<-data.frame(cbind(params_reg_unp50,d_unp50,interpret_unp50))


# 2.3 Overall level of avoidance

# Setting up data frame for GLM
df_unpcons<-dff4b3[,c(2,3,4,11,13,14,15)]

reg_unpcons_all<-glm(cbind(unpalatable_consumed,palatable_consumed)~.^2,data = df_unpcons,family = "binomial")

# Model selection
selection<-(step(reg_unpcons_all,scope = .~.^3))
summary(selection)

# Best model
reg_unpcons_best<-glm(formula = cbind(unpalatable_consumed, palatable_consumed) ~ 
                        sex + group + exp + pc1 + sex:group + sex:exp + sex:pc1 + 
                        exp:pc1 + sex:exp:pc1, family = "binomial", data = df_unpcons)

# Checking for overdispersion
overdisp_fun(reg_unpcons_best)

# chisq      ratio        rdf          p 
# 47.8444985  0.7249166 66.0000000  0.9549152 

# McFadden's pseudo R^2
with(summary(reg_unpcons_best), 1 - deviance/null.deviance)

# Model validation
plot(E<-simulateResiduals(reg_unpcons_best))

par(mfrow=c(2,3))
plotResiduals(E,as.factor(df_unpcons$sex))
plotResiduals(E,as.factor(df_unpcons$group))
plotResiduals(E,as.factor(df_unpcons$exp))
plotResiduals(E,as.factor(df_unpcons$population))
plotResiduals(E,df_unpcons$pc1)

# Model summary
summary(reg_unpcons_best)

# Effect size
params_reg_unpcons<-parameters::model_parameters(reg_unpcons_best) 
d_unpcons<-oddsratio_to_d(params_reg_unpcons$Coefficient,log = TRUE)
interpret_unpcons<-interpret_cohens_d(d_unpcons,rules = "gignac2016")
effectsize_unpcons<-data.frame(cbind(params_reg_unpcons,d_unpcons,interpret_unpcons))


# 2.4 Foraging latency

# Setting up the data frame for zero-inflated gamma regression
df.zeroif<-dff4b3[,c(1,2,3,4,6,14,15)]
outliers::grubbs.test(df.zeroif$latency)

# Excluding 9 latency outliers
exclude<-df.tweedie[df.tweedie$latency %in% c(28,33,39,45,46,54,62,65,100),]$id
df.zeroif2<-subset(df.zeroif,subset = !(df.zeroif$id %in% exclude))

# zero-inflated gamma regression
reg_latency_full<-glmmTMB(latency~(sex+group+exp+population)*pc1,ziformula = ~.,family = ziGamma(link = "log"),data = df.zeroif2)

# Model selection
selection<-step(reg_latency_full,direction = "both")
summary(selection)

# Best model
reg_latency_best<-glmmTMB(latency ~ sex + population + pc1 + sex:pc1,data = df.tweedie2,ziformula = ~.,family = ziGamma(link = "log"))

# model validation
plot(E<-simulateResiduals(reg_latency_best))
plotResiduals(E,as.factor(df.zeroif2$sex))
plotResiduals(E,df.zeroif2$pc1)

# Model summary
summary(reg_latency_best)