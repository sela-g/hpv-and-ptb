# title: "HPV and PTB"
# author: "Arianne Albert, PhD edited for dataset (8-MAY-2023) by Sela Grays"
# date: "05/27/2020"
# output:
  # html_document

rm(list = ls())
getwd()
setwd("/Users/selagrays/dev/hpv-and-ptb")


knitr::opts_chunk$set(echo = FALSE, include = FALSE)
options(digits = 2)

library(Gmisc)
library(broman)
library(pander)
library(knitr)
library(dplyr)
library(MASS)
library(effects)
library(emmeans)
# library(rcompanion)
library(DescTools)
library(tidyr)
library(sjPlot)
library(here)
library(readxl)
library(ggpubr)
library(blm)

theme_set(theme_pubr)

theme_pubr2 <- theme_pubr() +
  theme(panel.grid.major = element_line(linetype = 2, color = "grey"),
        text = element_text(size = 14))

theme_set(theme_pubr2)

cr2 <- function (formula, data, weights = NULL, na.action = na.omit) 
{
  risk <- function(start, stop, y, group, weights) {
    lower.upper <- quantile(group, c(start, stop))
    index <- which(group >= lower.upper[1] & group <= lower.upper[2])
    sum(y[index])/sum(weights[index])
  }
  covariate <- function(start, stop, group, weights) {
    lower.upper <- quantile(group, c(start, stop))
    index <- which(group >= lower.upper[1] & group <= lower.upper[2])
    sum(group[index] * weights[index])/sum(weights[index])
  }
  data <- data[, all.vars(formula)]
  print(dim(data))
  if (is.null(weights)) 
    weights <- rep(1, dim(data)[1])
  data$weights <- weights
  data <- na.action(data)
  y <- data[, all.vars(formula)[1]]
  print(y)
  group <- data[, all.vars(formula)[2]]
  Y <- mapply(risk, start = seq(0, 0.8, by = 0.01), stop = seq(0.2, 
                                                               1, by = 0.01), MoreArgs = list(group = group, y = y, 
                                                                                              weights = data$weights))
  X <- mapply(covariate, start = seq(0, 0.8, by = 0.01), stop = seq(0.2, 
                                                                    1, by = 0.01), MoreArgs = list(group = group, weights = data$weights))
  data.frame(risk = Y, x = X)
}

# ## vaccine data
vacc <- read.csv(here("perinatal_panorama_joined_external.csv"))
head(vacc)

## are there multiple rows per mom (to account for multiples)
any(duplicated(vacc$mother_study_id))
# TRUE


## also for multiple vaccine doses? (missing in dataset 16 Jun 2023)
# any(duplicated(vacc$BCCDC_ID))
# FALSE


## link to other data
## labour and delivery data
del.dat <- read.csv(here("R2018007_Multiple_Labours.csv"))
head(del.dat)

## mom data
mom.dat <- read.csv(here("R2018007_Mother.csv"))
head(mom.dat)

## baby data
baby.dat <- read.csv(here("R2018007_Baby.csv"))
head(baby.dat)
# not sure there's much of use here. The GA is more accurate in the delivery data


#del_imms <- read.csv(here("DeID Imms Records.csv")) # not sure that these match anymore

# add placeholders for missing vaccine dates
vacc[vacc$Immunization.Date == "", "Immunization.Date"] = "1000-01-01"

# manually calculate dose counts, subtract doses that occurred after giving birth
vacc <- vacc %>% group_by(mother_study_id) %>% 
  mutate(dose_counts = case_when(!is.na(Age.at.dose) ~ n_distinct(Immunization.Date) - sum(
    # vaccinated year(s) after giving birth
    strftime(Immunization.Date, format = "%Y") > baby_delivered_year | 
    # vaccinated month(s) after giving birth 
    # (NOTE: NO WAY OF KNOWING IF VACCINATED WITHIN SAME MONTH AS BIRTH, NO OCCURENCES IN DATASET (14 JUL 2023))
    (strftime(Immunization.Date, format = "%Y") > baby_delivered_year & strftime(Immunization.Date, format = "%m") > baby_delivered_month)), TRUE ~ 0))

# vaccination status is valid if there is vaccination data and all doses occurred after giving birth
vacc <- vacc %>%
  mutate(HPV.Vaccination.Status = case_when(
    !is.na(Age.at.dose) & dose_counts > 0 ~ "Valid",
    TRUE ~ "NA"
  ))
levels(vacc$HPV.Vaccination.Status) <- c("Valid", "NA")

# take the last row
vacc_mom <- vacc %>% slice(1)
# 5399 women

describeFactors(vacc_mom$m_num_births)
# 1 "5,786 (99.0%)"
# 2 "56 (1.0%)" 
# 28 sets of twins

# total number of babies
vacc %>% group_by(baby_study_id) %>% count %>% nrow()

describeFactors(vacc_mom$HPV.Vaccination.Status)

describeFactors(vacc_mom$dose_counts)
# Invalid "4 (0.1%)"     
# NA      "5,030 (93.2%)"
# Valid   "276 (5.1%)"   
# Missing "89 (1.6%)" 

vacc_mom <- vacc_mom %>% 
  mutate(vaccine = case_when(
    HPV.Vaccination.Status == "NA" ~ "No",
    HPV.Vaccination.Status == "Valid" ~ "Yes"
  ))

# the following columns are missing as of 08 Mar 2023
#vacc_mom$BCCDC_ID <- as.character(vacc_mom$BCCDC_ID)
#del_imms$BCCDC_Imms_Study_ID <- as.character(del_imms$BCCDC_Imms_Study_ID)

#vacc_test <- left_join(del_imms, vacc_mom, by = c("BCCDC_Imms_Study_ID" = "BCCDC_ID"))

# vacc_mom here should be vacc_test if/when the BCCDC_ID column is restored
vacc_test_one <- vacc_mom %>% 
  group_by(mother_study_id) %>%
  slice(1)

vacc_test_one <- vacc_test_one %>% 
  mutate(vacc_status = case_when(
    HPV.Vaccination.Status == "Valid" ~ "Yes",
    TRUE ~ "No"
  ))

describeFactors(vacc_test_one$vacc_status)

# del_imms is not present in combined dataset
# describeFactors(del_imms$HPV.Dose.Number)


#head(vacc2)
#head(del.dat) # del.dat does not have the linking id
#head(mom.dat) # this one links Imms study ID with Mother study ID

#dim(mom.dat) #[1] 5399   24
#dim(del.dat) #[1] 5447   15
# so there is more delivery data than mom data likely because some moms have more than one delivery?

# mom.del.merge <- merge(mom.dat, del.dat, all.y = TRUE, by.x = "BCCDC_Mother_Study_ID", by.y = "BCCDC_mother_study_ID")
# head(mom.del.merge)
# 
# # now can be merged with the HPV vaccine data
# vacc_mom <- merge(vacc2, mom.del.merge, all.y = TRUE, by.x = "BCCDC_Imms_Study_ID", by.y = "BCCDC_Imms_Study_ID")
# head(vacc_mom)
# 
# vacc_mom2 <- merge(vacc_mom, vacc_add2, all.x = TRUE, by.x = "BCCDC_Imms_Study_ID", by.y = "BCCDC_ID")
# head(vacc_mom)
# 
# # remove all the ones with IMMS ID == 0?? No I think we keep these as vaccine unknown status
# # vacc_mom <- vacc_mom[-which(vacc_mom$BCCDC_Imms_Study_ID == 0), ]
# # 5,337 deliveries
# 
# length(levels(factor(vacc_mom$BCCDC_Imms_Study_ID)))
# # 5,072 women
# 
# # remove anyone whose HPV dose is NA because we can't know their vaccination status?
# # AGAIN NO BECAUSE ITS JUST THAT WE DONT KNW THEIR STATUS AND WE SHOULD KEEP THE IN THE CONTROL/COMPARISON GROUP
# 
# # vacc_mom <- vacc_mom[-which(is.na(vacc_mom$HPV.Dose.Number)), ]
# # this leaves 4,281 deliveries because I think some of the ones in the vaccine file couldn't be matched to the delivery data?
# ## there is no study ID 2.... or 3 or 5??
# # this step makes sense
# 
# length(levels(factor(vacc_mom$BCCDC_Mother_Study_ID)))
# # 5399 women
# 
# # not sure how this number is larger than what was in the vaccine file?
# length(levels(factor(vacc_mom$BCCDC_Imms_Study_ID)))
# # 5073
# 
# # make new variable that is vaccinated vs not
# 
# vacc_mom$HPV.Vaccination.Status <- factor(ifelse(vacc_mom$HPV.Dose.Number == 0 | is.na(vacc_mom$HPV.Dose.Number), "no/unknown", "yes"))

# 
# ## lets also make a vaccine dose variable
# vacc_mom$vaccine.dose <- factor(case_when(
#   vacc_mom$HPV.Dose.Number == 0 ~ "none/unknown",
#   is.na(vacc_mom$HPV.Dose.Number) ~ "none/unknown", 
#   vacc_mom$HPV.Dose.Number == 1 ~ "one",
#   vacc_mom$HPV.Dose.Number == 2 ~ "two", 
#   TRUE ~"3 or more"), levels = c("none/unknown", "one", "two", "3 or more"))
# 

describeFactors(vacc_mom$dose_counts)
# # none/unknown "2,803 (51.5%)"
# # one          "198 (3.6%)"   
# # two          "451 (8.3%)"   
# # 3 or more    "1,995 (36.6%)"


getDescriptionStatsBy(vacc_mom$final_ga, vacc_mom$HPV.Vaccination.Status)
#           no/unknown           yes                 
# Mean (SD) "38.4 (&plusmn;2.2)" "38.3 (&plusmn;2.1)"
# Missing   "4 (0.1%)"           "3 (0.1%)" 

## make a preterm variable? ####
# what is the range of GA?
summary(vacc_mom$final_ga[which(vacc_mom$final_ga < 37)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 20.0    34.0    35.0    34.2    36.0    36.0

vacc_baby <- vacc %>% group_by(baby_study_id) %>% 
  mutate(PTB = factor(ifelse(final_ga < 37, "preterm", "term"), levels = c("term", "preterm"))) %>%
  slice(1)

vacc_mom <- vacc_mom %>%
  mutate(PTB = factor(ifelse(final_ga < 37, "preterm", "term"), levels = c("term", "preterm")))

describeFactors(vacc_baby$dose_counts)

getDescriptionStatsBy(vacc_baby$final_ga, vacc_baby$HPV.Vaccination.Status)

summary(vacc_mom$final_ga[which(vacc_mom$final_ga < 37)])

summary(vacc_mom$final_ga)

describeFactors(vacc_baby$HPV.Vaccination.Status)
# no/unknown "2,803 (51.5%)"
# yes        "2,644 (48.5%)"

describeFactors(vacc_baby$PTB)
# term    "4,794 (88.0%)"
# preterm "646 (11.9%)"  
# Missing "7 (0.1%)"

getDescriptionStatsBy(vacc_baby$PTB, vacc_baby$HPV.Vaccination.Status)
#         no/unknown      yes            
# term    "2,482 (88.5%)" "2,312 (87.4%)"
# preterm "317 (11.3%)"   "329 (12.4%)"  
# Missing "4 (0.1%)"      "3 (0.1%)" 

vacc_baby <- vacc_baby %>% 
  mutate(dose_counts = case_when(dose_counts < 1 | !is.numeric(dose_counts) ~ "none/unknown",
                                 dose_counts == 1 ~ "one", 
                                 dose_counts == 2 ~ "two", 
                                 dose_counts >= 3 ~ "3 or more"))
levels(vacc_baby$dose_counts) <- c("none/unknown", "one", "two", "3 or more")

getDescriptionStatsBy(vacc_baby$PTB, vacc_baby$dose_counts)
#         none/unknown    one           two           3 or more      
# term    "2,482 (88.5%)" "177 (89.4%)" "384 (85.1%)" "1,751 (87.8%)"
# preterm "317 (11.3%)"   "21 (10.6%)"  "66 (14.6%)"  "242 (12.1%)"  
# Missing "4 (0.1%)"      "0 (0.0%)"    "1 (0.2%)"    "2 (0.1%)"
# 08-Mar-2023
#         none/unknown    one           two           3 or more
# term    "2,158 (88.7%)" "206 (88.8%)" "437 (85.4%)" "1,993 (87.8%)"
# preterm "273 (11.2%)"   "26 (11.2%)"  "74 (14.5%)"  "273 (12.0%)"
# Missing "3 (0.1%)"      "0  (0.0%)"   "1 (0.2%)"    "3 (0.1%)"


### will want to control for things like age, and previous PTB, etc?

# split into different ranges
vacc_baby$ptb.3cat <- factor(ifelse(vacc_baby$final_ga < 32, "20-32", 
                                    ifelse(vacc_baby$final_ga < 35, "32-35", 
                                           ifelse(vacc_baby$final_ga < 37, "35-37", "term"))))
vacc_mom$ptb.3cat <- factor(ifelse(vacc_mom$final_ga < 32, "20-32", 
                                    ifelse(vacc_baby$final_ga < 35, "32-35", 
                                           ifelse(vacc_baby$final_ga < 37, "35-37", "term"))))
getDescriptionStatsBy(vacc_baby$ptb.3cat, vacc_baby$HPV.Vaccination.Status)
#         no/unknown      yes            
# 20-32   "44 (1.6%)"     "33 (1.2%)"    
# 32-35   "82 (2.9%)"     "71 (2.7%)"    
# 35-37   "191 (6.8%)"    "225 (8.5%)"   
# term    "2,482 (88.5%)" "2,312 (87.4%)"
# Missing "4 (0.1%)"      "3 (0.1%)"  


getDescriptionStatsBy(factor(vacc_baby$premature), vacc_baby$HPV.Vaccination.Status)
# these are previous preterm deliveries? I think so
#   no/unknown      yes            
# 0 "2,736 (97.6%)" "2,571 (97.2%)"
# 1 "59 (2.1%)"     "70 (2.6%)"    
# 2 "8 (0.3%)"      "2 (0.1%)"     
# 3 "0 (0.0%)"      "1 (0.0%)"


# make a category
vacc_mom$prev.ptb <- factor(ifelse(vacc_mom$premature == 0, "None", "At least one"), levels = c("None", "At least one"))
vacc_baby$prev.ptb <- factor(ifelse(vacc_baby$premature == 0, "None", "At least one"), levels = c("None", "At least one"))
describeFactors(vacc_baby$prev.ptb)
# None         "5,307 (97.4%)"
# At least one "140 (2.6%)"


getDescriptionStatsBy(factor(vacc_baby$r_substance_use), vacc_baby$HPV.Vaccination.Status)
#         no/unknown      yes            
# 1       "426 (15.2%)"   "563 (21.3%)"  
# Missing "2,377 (84.8%)" "2,081 (78.7%)"
# the missings are no I think


## thinking about other substance use definition
vacc_baby <- vacc_baby %>% 
  mutate(subs_2 = case_when(
    r_heroin == 1 | r_cocaine == 1 | r_methadone == 1 | r_solvents == 1 | r_rx == 1 | r_marijuana == 1 | r_other_drug == 1 | r_unk_drug == 1 ~ 1,
    TRUE ~ 0
  ))

describeFactors(vacc_baby$subs_2)
describeFactors(vacc_baby$r_marijuana)
describeFactors(vacc_baby$r_heroin)
describeFactors(vacc_baby$r_cocaine)
describeFactors(vacc_baby$r_methadone)
describeFactors(vacc_baby$r_solvents)
describeFactors(vacc_baby$r_rx)
describeFactors(vacc_baby$r_other_drug)
describeFactors(vacc_baby$r_unk_drug)

vacc_baby$subs_use <- replace(vacc_baby$r_substance_use, which(is.na(vacc_baby$r_substance_use)), 0)
vacc_baby$subs_use <- factor(vacc_baby$subs_use)
vacc_mom$subs_use <- replace(vacc_mom$r_substance_use, which(is.na(vacc_mom$r_substance_use)), 0)
vacc_mom$subs_use <- factor(vacc_mom$subs_use)
levels(vacc_baby$subs_use)


getDescriptionStatsBy(factor(vacc_baby$r_heroin), vacc_baby$HPV.Vaccination.Status)
#         no/unknown      yes            
# 1       "31 (1.1%)"     "28 (1.1%)"    
# Missing "2,772 (98.9%)" "2,616 (98.9%)"


getDescriptionStatsBy(factor(vacc_baby$r_cocaine), vacc_baby$HPV.Vaccination.Status)
#         no/unknown      yes            
# 1       "57 (2.0%)"     "65 (2.5%)"    
# Missing "2,746 (98.0%)" "2,579 (97.5%)"


getDescriptionStatsBy(factor(vacc_baby$r_methadone), vacc_baby$HPV.Vaccination.Status)
#         no/unknown      yes            
# 1       "27 (1.0%)"     "25 (0.9%)"    
# Missing "2,776 (99.0%)" "2,619 (99.1%)"


getDescriptionStatsBy(factor(vacc_baby$r_rx), vacc_baby$HPV.Vaccination.Status)
#         no/unknown      yes            
# 1       "11 (0.4%)"     "26 (1.0%)"    
# Missing "2,792 (99.6%)" "2,618 (99.0%)"


# getDescriptionStatsBy(vacc_baby$r_alc_flg, vacc_baby$HPV.Vaccination.Status)
# # not sure about this one. Will need to check data dictionary
# 
# getDescriptionStatsBy(factor(vacc_baby$second_hand_smoke), vacc_baby$HPV.Vaccination.Status)
# # not sure about this one. Will need to check data dictionary. Don't use has been removed from PSBC for low completion
# 
# getDescriptionStatsBy(factor(vacc_baby$smoker_type_cd), vacc_baby$HPV.Vaccination.Status)
# # not sure about this one. Will need to check data dictionary

vacc_baby$smoke.cat <- vacc_baby$smoker_type_cd
# NULL is a no
levels(vacc_baby$smoke.cat)[1] <- "N"
vacc_baby <- vacc_baby %>% rowwise() %>% mutate(smoke_flag = ifelse(smoke.cat == "C", 1, 0))
vacc_baby$smoke.cat <- factor(vacc_baby$smoke.cat, levels = c("N", "F", "C"), labels = c("Nonsmoker", "Quit prior to pregnancy", "Smoked during pregnancy"))
getDescriptionStatsBy(factor(vacc_baby$smoke.cat), vacc_baby$HPV.Vaccination.Status)


vacc_mom$smoke.cat <- vacc_mom$smoker_type_cd
# NULL is a no
levels(vacc_mom$smoke.cat)[1] <- "N"
vacc_mom <- vacc_mom %>% rowwise() %>% mutate(smoke_flag = ifelse(smoke.cat == "C", 1, 0))
vacc_mom$smoke.cat <- factor(vacc_mom$smoke.cat, levels = c("N", "F", "C"), labels = c("Nonsmoker", "Quit prior to pregnancy", "Smoked during pregnancy"))
getDescriptionStatsBy(factor(vacc_mom$smoke.cat), vacc_mom$HPV.Vaccination.Status)
# 08-MAR-2023
#NA	Valid
#Nonsmoker	650 (27.0%)	579 (19.4%)
#Quit prior to pregnancy	334 (13.9%)	570 (19.1%)
#moked during pregnancy	480 (19.9%)	749 (25.1%)
#Missing	946 (39.3%)	1,091 (36.5%)

options(digits = 10)
pre.1 <- glm(PTB ~ HPV.Vaccination.Status, data = vacc_baby, family = "binomial")
summary(pre.1)
# Coefficients:
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -2.05792    0.05964 -34.504   <2e-16 ***
# HPV.Vaccination.Statusyes  0.10811    0.08384   1.289    0.197 


plot(Effect(pre.1, focal.predictors = c("HPV.Vaccination.Status")))

### with previous preterm birth
pre.2 <- glm(PTB ~ HPV.Vaccination.Status + prev.ptb, data = vacc_baby, family = "binomial")
summary(pre.2)
# Coefficients:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -2.11761    0.06089 -34.780   <2e-16 ***
# HPV.Vaccination.Statusyes     0.10074    0.08449   1.192    0.233    
# prev.ptbAt least one  1.50793    0.18100   8.331   <2e-16 ***



plot(Effect(pre.2, focal.predictors = c("prev.ptb")))
plot(Effect(pre.2, focal.predictors = c("HPV.Vaccination.Status", "prev.ptb")))

# subs_use
pre.3 <- glm(PTB ~ HPV.Vaccination.Status + prev.ptb + subs_use, data = vacc_baby, family = "binomial")
summary(pre.3)
# Coefficients:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -2.18839    0.06409 -34.147  < 2e-16 ***
# HPV.Vaccination.Statusyes     0.07279    0.08495   0.857    0.392    
# prev.ptbAt least one  1.52640    0.18161   8.405  < 2e-16 ***
# subs_use1             0.40924    0.10113   4.047 5.19e-05 ***



plot(Effect(pre.3, focal.predictors = c("subs_use")))
plot(Effect(pre.3, focal.predictors = c("subs_use", "prev.ptb")))

exp(cbind(pre.3$coefficients, confint(pre.3)))
#                                     2.5 %    97.5 %
# (Intercept)          0.1120973 0.09867723 0.1268702
# HPV.Vaccination.Statusyes    1.0755000 0.91051374 1.2704542
# prev.ptbAt least one 4.6015776 3.20488942 6.5416378
# subs_use1            1.5056749 1.23196183 1.8317455


## baby year of birth?

# # smoking
# vacc_baby$smoker_type_cd
# vacc_baby$smoke.cat <- factor(vacc_baby$smoke.cat, levels = c())

pre.4 <- glm(PTB ~ HPV.Vaccination.Status + prev.ptb + smoke_flag, data = vacc_baby, family = "binomial")
summary(pre.4)

drop1(pre.4, test = "Chi")
#                Df Deviance    AIC    LRT  Pr(>Chi)    
# <none>              3902.6 3912.6                     
# HPV.Vaccination.Status  1   3903.8 3911.8  1.205    0.2723    
# prev.ptb        1   3960.4 3968.4 57.802 2.899e-14 ***
# smoke.cat       2   3905.2 3911.2  2.545    0.2801  


pre.all <- glm(PTB ~ HPV.Vaccination.Status + prev.ptb + subs_use + smoke.cat, data = vacc_baby, family = "binomial")
summary(pre.all)
#08-MAR@-2023:
#                                 Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                       -2.2577     0.1093  -20.65  < 2e-16 ***
#  HPV.Vaccination.StatusValid        0.1198     0.1126    1.06     0.29    
#  prev.ptbAt least one               1.3181     0.2489    5.29  1.2e-07 ***
#  subs_use1                          0.4988     0.1270    3.93  8.6e-05 ***
#  smoke.catQuit prior to pregnancy  -0.1657     0.1473   -1.13     0.26    
#  smoke.catSmoked during pregnancy  -0.0227     0.1360   -0.17     0.87    


# using blm package to get a sense of if the odds ratios and logistic regression are properly estimating the risks
vacc_baby$PTB2 <- ifelse(vacc_baby$PTB == "preterm", 1, 0)

fit1 <- blm(PTB2 ~ HPV.Vaccination.Status, data = vacc_baby)
summary(fit1)

coef(fit1)*100

confint(fit1)*100
#                        Est.      Lower     Upper
# (Intercept)       11.325521 10.3383405 12.312701
# HPV.Vaccination.Statusyes  1.131907 -0.3405433  2.604357


# TODO: Is PTB.spont supposed to be preterm *and* spontaneous? then why denote "term"
# and not 0 1 etc?... as of 05Jul2023, I still don't genuinely understand what.. is happening here. 
levels(vacc_baby$PTB)
levels(as.factor(vacc_baby$m_labour_type))
vacc_baby <- vacc_baby %>% rowwise() %>% 
  mutate(PTB.spont = if_else(PTB == "term", 0,
                                           case_match(m_labour_type, "spontaneous" ~ 1,
                                                                    c("induced", "none") ~ 0)),
         PTB.spont2 = if_else(PTB == "term", "Term",
                              case_match(m_labour_type, "spontaneous" ~ "Spontaneous preterm",
                                         c("induced", "none") ~ "Iatrogenic preterm")))
levels(vacc_baby$PTB.spont2) <- c("Spontaneous preterm", "Iatrogenic preterm", "Term")

vacc_mom <- vacc_mom %>% rowwise() %>% 
  mutate(PTB.spont = if_else(PTB == "term", 0,
                             case_match(m_labour_type, "spontaneous" ~ 1,
                                        c("induced", "none") ~ 0)),
         PTB.spont2 = if_else(PTB == "term", "Term",
                              case_match(m_labour_type, "spontaneous" ~ "Spontaneous preterm",
                                         c("induced", "none") ~ "Iatrogenic preterm")))
levels(vacc_mom$PTB.spont2) <- c("Spontaneous preterm", "Iatrogenic preterm", "Term")

spon.1 <- glm(PTB.spont ~ HPV.Vaccination.Status, data = vacc_baby, family = "binomial")
summary(spon.1)
# Coefficients:
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -2.31149    0.06684 -34.581   <2e-16 ***
# HPV.Vaccination.Statusyes  0.02094    0.09578   0.219    0.827 

plot(Effect(spon.1, focal.predictors = c("HPV.Vaccination.Status")))

### with previous preterm birth
spon.2 <- glm(PTB.spont ~ HPV.Vaccination.Status + prev.ptb, data = vacc_baby, family = "binomial")
summary(spon.2)
# Coefficients:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -2.37611    0.06837 -34.756  < 2e-16 ***
# HPV.Vaccination.Statusyes     0.01346    0.09650   0.140    0.889    
# prev.ptbAt least one  1.56946    0.19680   7.975 1.52e-15 ***

plot(Effect(spon.2, focal.predictors = c("prev.ptb")))
plot(Effect(spon.2, focal.predictors = c("HPV.Vaccination.Status", "prev.ptb")))

# subs_use
spon.3 <- glm(PTB.spont ~ HPV.Vaccination.Status + prev.ptb + subs_use, data = vacc_baby, family = "binomial")
summary(spon.3)
# Coefficients:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -2.43919    0.07199 -33.883  < 2e-16 ***
# HPV.Vaccination.Statusyes    -0.01111    0.09696  -0.115   0.9088    
# prev.ptbAt least one  1.58768    0.19733   8.046 8.56e-16 ***
# subs_use1             0.36967    0.11645   3.175   0.0015 **  

plot(Effect(spon.3, focal.predictors = c("subs_use")))
plot(Effect(spon.3, focal.predictors = c("subs_use", "prev.ptb")))

exp(cbind(spon.3$coefficients, confint(spon.3)))
#                                      2.5 %    97.5 %
# (Intercept)          0.08723105 0.07556157 0.1002079
# HPV.Vaccination.Statusyes    0.98894972 0.81757463 1.1958651
# prev.ptbAt least one 4.89239492 3.29278441 7.1517987
# subs_use1            1.44725830 1.14767960 1.8123377

### year of birth
summary(vacc_baby$baby_delivered_year)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2015    2015    2016    2016    2016    2017 

vacc_baby$baby_delivered_year <- factor(vacc_baby$baby_delivered_year)
vacc_baby$PTB.spont <- factor(vacc_baby$PTB.spont)
getDescriptionStatsBy(vacc_baby$PTB.spont, vacc_baby$baby_delivered_year)
#           2015            2016            2017         
# term      "1,765 (90.9%)" "2,287 (90.6%)" "742 (90.8%)"
# spont.ptb "174 (9.0%)"    "231 (9.2%)"    "75 (9.2%)"  
# Missing   "2 (0.1%)"      "5 (0.2%)"      "0 (0.0%)" 


# year
spon.4 <- glm(PTB.spont ~ HPV.Vaccination.Status + prev.ptb + subs_use + baby_delivered_year, data = vacc_baby, family = "binomial")
summary(spon.4)
# Coefficients:
#                          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             -2.441596   0.096793 -25.225  < 2e-16 ***
# HPV.Vaccination.Statusyes       -0.010801   0.096998  -0.111  0.91134    
# prev.ptbAt least one     1.587712   0.197553   8.037 9.22e-16 ***
# subs_use1                0.370231   0.116509   3.178  0.00148 ** 
# baby_delivered_year2016 -0.001794   0.106231  -0.017  0.98653    
# baby_delivered_year2017  0.019252   0.146219   0.132  0.89525 


getDescriptionStatsBy(vacc_baby$ptb.3cat, vacc_baby$HPV.Vaccination.Status, statistics = TRUE)
#         no/unknown      yes             P-value
# 20-32   "33 (1.2%)"     "24 (0.9%)"     "0.31" 
# 32-35   "65 (2.4%)"     "50 (2.0%)"     ""     
# 35-37   "148 (5.4%)"    "160 (6.3%)"    ""     
# term    "2,482 (90.8%)" "2,312 (90.7%)" ""     
# Missing "4 (0.1%)"      "3 (0.1%)"      "" 

CochranArmitageTest(table(vacc_baby$HPV.Vaccination.Status, vacc_baby$ptb.3cat))
# p-value = 0.5507
#08-Mar-2023:
#p-value = 0.5447

getDescriptionStatsBy(vacc_baby$ptb.3cat, vacc_baby$dose_counts, statistics = TRUE)
#         none/unknown    one           two           3 or more       P-value
# 20-32   "33 (1.2%)"     "0 (0.0%)"    "2 (0.5%)"    "22 (1.1%)"     "0.28" 
# 32-35   "65 (2.4%)"     "4 (2.1%)"    "11 (2.5%)"   "35 (1.8%)"     ""     
# 35-37   "148 (5.4%)"    "9 (4.7%)"    "35 (8.1%)"   "116 (6.0%)"    ""     
# term    "2,482 (90.8%)" "177 (93.2%)" "384 (88.7%)" "1,751 (90.9%)" ""     
# Missing "4 (0.1%)"      "0 (0.0%)"    "1 (0.2%)"    "2 (0.1%)"      ""  


spon.all <- glm(PTB.spont ~ HPV.Vaccination.Status + prev.ptb + subs_use + smoke.cat, data = vacc_baby, family = "binomial")
summary(spon.all)

vacc_mom$gravida
describeMedian(vacc_mom$gravida)

vacc_baby$grav.cat <- factor(case_when(
  vacc_baby$gravida == 1 ~ "One",
  vacc_baby$gravida <= 3 ~ "Two to three",
  vacc_baby$gravida >3 ~ "Four or more"
),
levels = c("One", "Two to three", "Four or more"))
vacc_mom$grav.cat <- factor(case_when(
  vacc_mom$gravida == 1 ~ "One",
  vacc_mom$gravida <= 3 ~ "Two to three",
  vacc_mom$gravida >3 ~ "Four or more"
),
levels = c("One", "Two to three", "Four or more"))
describeFactors(vacc_baby$grav.cat)

vacc_baby$prev.ptb

vacc_baby$m_bmi_no <- replace(vacc_baby$m_bmi_no, vacc_baby$m_bmi_no == 9999.99, NA)

vacc_baby$subs_use
describeFactors(vacc_baby$subs_use)

vacc_baby$m_mode_del
describeFactors(vacc_baby$m_mode_del)

vacc_baby$final_ga
getDescriptionStatsBy(vacc_baby$final_ga, vacc_baby$HPV.Vaccination.Status)

vacc_baby$ptb.3cat

#vacc_mom$vaccine.dose


# vacc_mom[, varname] returns array, need to use vacc_mom[[varname]] to return vector
getT1Stat <- function(varname, digits=0, useNA = "ifany"){
  getDescriptionStatsBy(vacc_baby[[varname]], 
                        vacc_baby$HPV.Vaccination.Status, 
                        add_total_col=TRUE,
                        show_all_values=TRUE, 
                        hrzl_prop=FALSE,
                        statistics= FALSE, 
                        html=TRUE,
                        useNA = useNA,
                        digits=digits,
                        header_count = TRUE)
}

getT1Stat.median <- function(varname, digits=0, useNA = "ifany"){
  getDescriptionStatsBy(vacc_baby[[varname]], 
                        vacc_baby$HPV.Vaccination.Status, 
                        add_total_col=TRUE,
                        show_all_values=TRUE, 
                        hrzl_prop=FALSE,
                        statistics=FALSE, 
                        html=TRUE, 
                        useNA = useNA,
                        digits=digits,
                        header_count = TRUE,
                        continuous_fn = describeMedian)
}


table_data <- list()

# Get the basic stats
table_data[["Gravida"]] <- getT1Stat("grav.cat", 1)
table_data[["Previous preterm deliveries"]] <- getT1Stat.median("prev.ptb", 1)
table_data[["BMI"]] <- getT1Stat("m_bmi_no", 1)
table_data[["Smoking status"]] <- getT1Stat.median("smoke.cat", 1)
table_data[["Maternal substance use"]] <- getT1Stat("subs_use", 1)
table_data[["GA at delivery"]] <- getT1Stat.median("final_ga", 1)
table_data[["Preterm delivery"]] <- getT1Stat.median("ptb.3cat", 1)
table_data[["Spontaneous Preterm delivery"]] <- getT1Stat.median("PTB.spont2", 1)
table_data[["Mode of delivery"]] <- getT1Stat("m_mode_del", 1)
table_data[["Number of vaccine doses"]] <- getT1Stat("dose_counts", 1)



# Now merge everything into a matrix
# and create the rgroup & n.rgroup variabels
rgroup <- c()
n.rgroup <- c()
output_data <- NULL
for (varlabel in names(table_data)) {
  output_data <- rbind(output_data, table_data[[varlabel]])
  rgroup <- c(rgroup,
              varlabel)
  n.rgroup <- c(n.rgroup,
                nrow(table_data[[varlabel]]))
}

# Add a column spanner for the columns
cgroup <- c("", "Vaccination status")
n.cgroup <- c(1, 2)
colnames(output_data) <- gsub("[ ]*Vaccination status", "", colnames(output_data))

# htmlTable::htmlTable(output_data, align = "rrrr",
# rgroup = rgroup, n.rgroup = n.rgroup,
# css.rgroup.sep = "",
# cgroup = cgroup,
# n.cgroup = n.cgroup,
# rowlabel = "",
# caption = "Table 1. Demographic and clinical data summaries.",
# ctable = TRUE) %>% htmltools::html_print()

## Results

# Table 1 gives a summary of clinical and demographic variables.

htmlTable::htmlTable(output_data, align = "rrrr",
                     rgroup = rgroup, n.rgroup = n.rgroup,
                     css.rgroup.sep = "",
                     cgroup = cgroup,
                     n.cgroup = n.cgroup,
                     rowlabel = "",
                     caption = "Table 1. Demographic and clinical data summaries.",
                     ctable = TRUE)

tab_model(pre.1, pre.all, show.intercept = FALSE, show.r2 = FALSE, show.reflvl = FALSE, 
          terms = c("HPV.Vaccination.Status [NA, Valid]", "prev.ptb [None, At least one]", 
                    "subs_use [0, 1]", "smoke.cat [Smoked during pregnancy]"), 
          collapse.ci = FALSE, pred.labels = c("HPV vaccine", "Previous preterm delivery", "Substance use", "Smoked during pregnancy"), 
          dv.labels = c("Unadjusted model", "Adjusted model"), string.ci = "95% CI", string.est = "Risk Ratios")

pre.risk <- summary(Effect(pre.1, focal.predictors = "HPV.Vaccination.Status"))
pre.tab <- data.frame(risk = pre.risk$effect, vaccine = c("No/unknown", "Vaccinated"), low = pre.risk$lower, up = pre.risk$upper)

ggplot(pre.tab, aes(x = vaccine, y = risk)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.09, 0.14)) +
  geom_pointrange(aes(y = risk, ymin = low, ymax = up)) +
  xlab("HPV vaccine status") +
  ylab("Unadjusted absolute risk\nof preterm birth")

tab_model(spon.1, spon.all, show.intercept = FALSE, show.r2 = FALSE, show.reflvl = FALSE, 
          terms = c("HPV.Vaccination.Status [NA, Valid]", "prev.ptb [None, At least one]", 
                    "subs_use [0, 1]", "smoke.cat [Smoked during pregnancy]"),
          collapse.ci = FALSE, pred.labels = c("HPV vaccine", "Previous preterm delivery", "Substance use", "Smoked during pregnancy"), dv.labels = c("Unadjusted model", "Adjusted model"), string.ci = "95% CI", string.est = "Risk Ratios")

pre.risk2 <- summary(Effect(spon.1, focal.predictors = "HPV.Vaccination.Status"))
pre.tab2 <- data.frame(risk = pre.risk2$effect, vaccine = c("No/unknown", "Vaccinated"), low = pre.risk2$lower, up = pre.risk2$upper)

ggplot(pre.tab2, aes(x = vaccine, y = risk)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.055, 0.09)) +
  geom_pointrange(aes(y = risk, ymin = low, ymax = up)) +
  xlab("HPV vaccine status") +
  ylab("Unadjusted absolute risk\nof spontaneous preterm birth")


