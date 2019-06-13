### Final version of LME code accompanying the 2019 manuscript: 
### "Biomass encounter rates limit the size scaling of feeding interactions"
### By Daniel Barrios-O’Neill, Ruth Kelly & Mark C. Emmerson, 
###

### Model fitting by Daniel Barrios-O’Neill & Ruth Kelly


# Data is in .csv file - "benthic_functional_response_DB_withMAXMIN.csv"
# Metadata in word doc - "metadata_biomass_encounters_2019.docx"


### Packages required lme4, MuMIn, sjPlot for lme viz  
### install these before proceeding
### FYI: We used lme4 (version 1.1-21) and MuMIn for model selection (version 1.40.4)
### and sjPlot (version 2.4.1) in R version 3.3.3

#### Set-up ####

rm(list=ls())  ### clear workspace

#  setwd("C:/R/Barrios_ONeill_2019_biomass_e_rates") ## set your work directory

# load required packages

library("lme4")
library("MuMIn")
library("sjPlot")

## Read in dataset
d<-read.csv("benthic_functional_response_DB_withMAXMIN.csv",  na.strings = c("", "NA"))

#library(lme4) ## model fitting lme functions
#library(MuMIn) ## dredge, AICc etc

## Removing our categories of animals with insuffient data 
## SW: an *obligate* sit and wait predator (e.g. anemones)
## D: a deposit feeder, G = Grazers.
d2<-d[d$interaction!="D" & d$interaction!="SW" & d$interaction!="G",]

#### Combine predator type levels "Chondrichthyes" and "Osteichthyes" 
#### into one category named "fish"
levels(d2[,"pred_type"])[c(3,8)]<-c("fish") ### spotted a mistake here:

### the two fish types need to be combined, check:
d2<-droplevels(d2)
levels(d2[,"pred_type"])


##### Mixed model for Capture Rates ####

##### first dealing with slopes (capture rates)
global_cap<-lmer(log(slopes.g.m2.d)~((log(pred_g)+log(prey_g)+log(temp_C))^2
                +interaction)+ 
               (1|pred_type), ### random effect structure: pred type, i.e. higher level taxonomy
                data=d2, 
                REML = FALSE, 
                na.action = "na.fail")


### visualisation for model checking:
sjp.lmer(global_cap,"re", y.offset = .4) ## 
sjp.lmer(global_cap,"fe", y.offset = .4)
plot_model(global_cap,"diag") ## reasonable

### Running Null model with intercept only and comparing AICc value
### with the global model. Doing this to control for type 1 error rate after # Harrison et al (2018)
int_mod<-lmer(log(slopes.g.m2.d)~1+(1|pred_type), data=d2) ### intercept only
AICc(global_cap,int_mod) # grand.

#### Model selection
#### Here we use the function dredge to run all model subsets of the global model
#### and find "minimal adequate model" i.e. that with lowest AICc
options(na.action = "na.fail")
ms1 <- dredge(global_cap)


### Make neat AICc table for supp matt
ms1_table <- as.data.frame(ms1)
head(ms1_table)
names(ms1_table) <- c("Intercept", "Encounter strategy", "log(Predator mass)",        
                       "log(prey_mass)", "log(temperature)", "log(pred_mass)*log(prey_mass)",
                        "log(pred_mass):log(temperature)", "log(prey_mass):log(temperature)",
                        "df", "logLik", "AICc","delta", "weight") 

### remove unneeded columns
ms1_table <- ms1_table[,c(1:8,11)]

### write to .csv file
#write.csv(ms1_table, "AICc_capture_rates_13062019.csv", row.names = F)

#### Get top model by AICc and inspect 
bestmod<-(get.models(ms1, 1)[[1]])

summary(bestmod)
class(bestmod)
sjp.lmer(bestmod,"re", y.offset = .4) 
sjp.lmer(bestmod,"fe", y.offset = .4) ### 
plot_model(bestmod,"diag") ### acceptable


#### Fitting models for max feeding rates (MFR)

#### Prior to model fitting it is necessary to illiminate tiny values and 0's 
#### This affects 3 of 171 data points. 

#### eliminate 0s in handling times:
d3<-d2[d2$handling!="0",] ### eliminate 0 vals (n = 3)

### remove outlier large values in MFR
d3<-d3[d3$MFR<(exp(1))^3,] ## tiny vals (n = 20)


#####
globMFR<-lmer(log(MFR)~((log(pred_g)+log(prey_g)+log(temp_C))^2+interaction)+ 
               (1|pred_type), ### randon effect structure: pred type, i.e. higher level taxonomy
             data=d3, 
             REML = FALSE, 
             na.action = "na.fail")
### mod checks:
sjp.lmer(globMFR,"re", y.offset = .4) ## bigger REs here
sjp.lmer(globMFR,"fe", y.offset = .4) ## some basic body size effects are marginal.
plot_model(globMFR,"diag") ## not as good as slopes

int_MFR<-lmer(log(MFR)~1+(1|pred_type), data=d3) ### intercept only
AICc(globMFR,int_MFR) # grand. Doing this to control for type 1 error rate after # Harrison et al (2018)

#### dredge for min mod:
options(na.action = "na.fail")
ms2 <- dredge(globMFR)

bestMFR<-(get.models(ms2, 1)[[1]])

sjp.lmer(bestMFR,"re", y.offset = .4) ### REs still important here. 
sjp.lmer(bestMFR,"fe", y.offset = .4) ### makes sense
plot_model(bestMFR,"diag") ### 



### Make neat AICc table for supp matt
ms2_table <- as.data.frame(ms2)
head(ms2_table)
names(ms2_table) <- c("Intercept", "Encounter strategy", "log(Predator mass)",        
                      "log(prey_mass)", "log(temperature)", "log(pred_mass)*log(prey_mass)",
                      "log(pred_mass):log(temperature)", "log(prey_mass):log(temperature)",
                      "df", "logLik", "AICc","delta", "weight") 

### remove unneeded columns
ms2_table <- ms2_table[,c(1:8,11)]

### write to .csv file
write.csv(ms2_table, "AICc_MFR_13062019.csv", row.names = F)

