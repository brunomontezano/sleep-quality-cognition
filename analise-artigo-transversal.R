##########################################################
### Statistical analysis for article entitled          ###
### "Associations between general sleep quality        ###
### and measures of functioning and cognition in       ###
### subjects recently diagnosed with bipolar disorder" ###
### Author: Bruno Montezano                            ###
##########################################################

###############################
###### LIBRARY AND DATA LOADING
###############################

# Load libraries
library(tidyverse)
library(psych)
library(foreign)
library(FSA)
library(PMCMRplus)
#library(summarytools)
library(Hmisc)
library(ggpubr)

# Remove scientific notation
options(scipen = 35)

# Load database from file
setwd("/home/pepper/dox/uni/tcp2/")
df <- read.spss("banco-conversao-16-10-20.sav",
                to.data.frame=TRUE,
                quiet=TRUE,
                use.label.values=FALSE,
                stringAsFactors=FALSE)

# Creating a second dataframe with follow-up observations
fu <- df %>% filter(., !is.na(mora_t2))

# Creating variable, and calculating the WAIS Letter-Number subtest score
nl <- fu %>% 
      select(starts_with("numletr", ignore.case = TRUE)) %>% 
      mutate_all(~replace(., is.na(.), 0))

fu <- fu %>% 
      mutate(NUMLETR_total = rowSums(nl[1:21]))

# Article groups variable (BD, MDD in current episode, euthymic MDD)
fu <- fu %>% 
      mutate(grupos = case_when(
                           TB_erros == "conversão para TB" ~ 'Transtorno bipolar',
                           (miniA08AT_t2 == 1 |
                            miniA08ATPA_t2 == 1 |
                            miniA15b_t2 == 1)
                            &
                            TB_erros != "conversão para TB" ~ 'TDM em episódio atual',
                           (is.na(miniA08AT_t2) &
                            is.na(miniA08ATPA_t2) &
                            is.na(miniA15b_t2))
                            &
                            TB_erros != "conversão para TB" ~ 'TDM eutímico'
                      ))

# Import psychiatric/pharmacological treatment during life variable
df_tai <- read.spss("banco-conversao-t1-t2-mario2.sav",
                to.data.frame=TRUE,
                quiet=TRUE,
                use.label.values=FALSE,
                stringAsFactors=FALSE)

# Renaming REC variable name in Tai's database
colnames(df_tai)[colnames(df_tai) == "a02rec"] <- "rec"

# Creating subset with just REC, treatment, and marital status variables
df_tai_sub <- df_tai %>%
   select(rec, lifetimepsychiatricmed_t2, maritalstatusdic_t2)

# Joining subset that I just created to the follow-up dataframe
fu <- left_join(fu, df_tai_sub, by = "rec")

# Replace NA on suicide risk variable to no-risk subject
fu$riscodesuicidioatual <- fu$riscodesuicidioatual %>% 
                           replace_na("não")

# Create anxiety disorders variable

fu <- fu %>% 
      mutate(anxiety_dic = case_when(
                           miniF04c_t2 == 1 |
                           miniF04d_t2 == 1 |
                           miniP06_t2 == 1 |
                           miniG06_t2 == 1 |
                           miniG07_t2 == 1 ~ 1,
                           TRUE ~ 0)
                      )

# Creating BD subset
bd <- fu %>%
      filter(., grupos == 'BD')

# Creating MDD in current episode subset
mdd_cur <- fu %>% 
           filter(., grupos == 'MDD in current episode')

# Creating euthymic MDD subset
mdd_eut <- fu %>% 
           filter(., grupos == 'Euthymic MDD')

# scatterplots para apresentação do salao 2021
plot_fast <- fu %>%
  ggplot(aes(x = PSQI_total, y = FAST_total_t2, color = grupos)) +
  geom_jitter(alpha = 0.8) +
  labs(x = "Escore da PSQI", y = "Escore da FAST") +
  geom_smooth(method = "lm") +
  theme_light(base_size = 16) +
  scale_color_viridis_d(name = "Grupos", option = "plasma") +
  theme(legend.position = "none")

plot_cobra <- fu %>%
  ggplot(aes(x = PSQI_total, y = COBRA_soma_t2_16itens, color = grupos)) +
  geom_jitter(alpha = 0.8) +
  labs(x = "Escore da PSQI", y = "Escore da COBRA") +
  geom_smooth(method = "lm") +
  theme_light(base_size = 16) +
  scale_color_viridis_d(name = "Grupos", option = "plasma") +
  theme(legend.position = "none")

plot_numletr <- fu %>%
  ggplot(aes(x = PSQI_total, y = NUMLETR_total, color = grupos)) +
  geom_jitter(alpha = 0.8) +
  labs(x = "Escore da PSQI", y = "Escore do subteste da WAIS") +
  geom_smooth(method = "lm") +
  theme_light(base_size = 16) +
  scale_color_viridis_d(name = "Grupos", option = "plasma") +
  theme(legend.position = "none")

legend <- cowplot::get_legend(
  plot_fast +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_pronto <- cowplot::plot_grid(
  plot_fast, plot_cobra, plot_numletr,
  legend)

ggplot2::ggsave(filename = "~/tmp/plot_correlacoes.png",
  plot = plot_pronto,
  height = 3800,
  width = 5200,
  units = "px",
  dpi = 300)

# Test distribution

ggqqplot(bd$FAST_total_t2)
ggqqplot(bd$COBRA_soma_t2_16itens)
ggqqplot(bd$NUMLETR_total)

ggqqplot(mdd_eut$FAST_total_t2)
ggqqplot(mdd_eut$COBRA_soma_t2_16itens)
ggqqplot(mdd_eut$NUMLETR_total)

ggqqplot(mdd_cur$FAST_total_t2)
ggqqplot(mdd_cur$COBRA_soma_t2_16itens)
ggqqplot(mdd_cur$NUMLETR_total)

shapiro.test(fu$FAST_total_t2)
shapiro.test(fu$COBRA_soma_t2_16itens)
shapiro.test(fu$NUMLETR_total)

shapiro.test(bd$FAST_total_t2)
shapiro.test(bd$COBRA_soma_t2_16itens)
shapiro.test(bd$NUMLETR_total)

shapiro.test(mdd_eut$FAST_total_t2)
shapiro.test(mdd_eut$COBRA_soma_t2_16itens)
shapiro.test(mdd_eut$NUMLETR_total)

shapiro.test(mdd_cur$FAST_total_t2)
shapiro.test(mdd_cur$COBRA_soma_t2_16itens)
shapiro.test(mdd_cur$NUMLETR_total)
############################################################
###### TABLE 1 DATA: SOCIODEMOGRAPHIC AND CLINICAL VARIABLES
############################################################

### SEX ###

freq(bd$a03sexo_t2)
freq(mdd_cur$a03sexo_t2)
freq(mdd_eut$a03sexo_t2)
chisq.test(fu$a03sexo_t2, fu$grupos)

### AGE ###

idade <- fu %>% select(., a05idade_t2)
stby(data = idade,
     INDICES = fu$grupos, # Stratify by group FUN = descr, # Descriptive analysis
     stats = "common") # Most common stats

aov_idade <- aov(a05idade_t2 ~ grupos, data = fu)
summary(aov_idade)

### YEARS OF EDUCATION ###

escolaridade <- fu %>% select(., escolaridade)
stby(data = escolaridade,
     INDICES = fu$grupos, # Stratify by group
     FUN = descr, # Descriptive analysis
     stats = "common") # Most common stats

aov_escolaridade <- aov(escolaridade ~ grupos, data = fu)
summary(aov_escolaridade)

### STUDYING/WORKING ###

freq(bd$nemtrabnemestuda)
freq(mdd_cur$nemtrabnemestuda)
freq(mdd_eut$nemtrabnemestuda)
chisq.test(fu$nemtrabnemestuda, fu$grupos)

### MARITAL STATUS ###

freq(bd$maritalstatusdic_t2)
freq(mdd_cur$maritalstatusdic_t2)
freq(mdd_eut$maritalstatusdic_t2)
chisq.test(fu$maritalstatusdic_t2, fu$grupos)

### ETHNICITY ###

freq(bd$cordapele)
freq(mdd_cur$cordapele)
freq(mdd_eut$cordapele)
chisq.test(fu$cordapele, fu$grupos)

### SOCIOECONOMIC STATUS ###

freq(bd$abepdicotomica)
freq(mdd_cur$abepdicotomica)
freq(mdd_eut$abepdicotomica)
chisq.test(fu$abepdicotomica, fu$grupos)

### ALCOHOL ABUSE/DEPENDENCE ###

freq(bd$alcoolabudep)
freq(mdd_cur$alcoolabudep)
freq(mdd_eut$alcoolabudep)
chisq.test(fu$alcoolabudep, fu$grupos)

### ILLICIT SUBSTANCE ABUSE/DEPENDENCE ###

freq(bd$abudepoutrasdrogas)
freq(mdd_cur$abudepoutrasdrogas)
freq(mdd_eut$abudepoutrasdrogas)
chisq.test(fu$abudepoutrasdrogas, fu$grupos)

### CURRENT SUICIDE RISK ###

freq(bd$riscodesuicidioatual)
freq(mdd_cur$riscodesuicidioatual)
freq(mdd_eut$riscodesuicidioatual)
chisq.test(fu$riscodesuicidioatual, fu$grupos)

### LIFETIME PSYCHIATRIC MEDICATION ###

freq(bd$lifetimepsychiatricmed_t2)
freq(mdd_cur$lifetimepsychiatricmed_t2)
freq(mdd_eut$lifetimepsychiatricmed_t2)
chisq.test(fu$lifetimepsychiatricmed_t2, fu$grupos)

### HYPNOTICS/SEDATIVES ABUSE/DEPENDENCE ###

freq(bd$hipnoticosabudep)
freq(mdd_cur$hipnoticosabudep)
freq(mdd_eut$hipnoticosabudep)
chisq.test(fu$hipnoticosabudep, fu$grupos)

### PSQI ###

#quantile(bd$PSQI_total, probs = c(0.5,0.25,0.75))
#quantile(mdd_eut$PSQI_total, probs = c(0.5,0.25,0.75), na.rm = TRUE)
#quantile(mdd_cur$PSQI_total, probs = c(0.5,0.25,0.75), na.rm = TRUE)
#kruskal.test(PSQI_total ~ as.factor(grupos), data = fu)

summary(aov(PSQI_total ~ grupos, data = fu))

mean(bd$PSQI_total)
sd(bd$PSQI_total)
mean(mdd_eut$PSQI_total, na.rm = TRUE)
sd(mdd_eut$PSQI_total, na.rm = TRUE)
mean(mdd_cur$PSQI_total, na.rm = TRUE)
sd(mdd_cur$PSQI_total, na.rm = TRUE)

### FAST ###

#quantile(bd$FAST_total_t2, probs = c(0.5,0.25,0.75))
#quantile(mdd_eut$FAST_total_t2, probs = c(0.5,0.25,0.75), na.rm = TRUE)
#quantile(mdd_cur$FAST_total_t2, probs = c(0.5,0.25,0.75), na.rm = TRUE)
#kruskal.test(FAST_total_t2 ~ as.factor(grupos), data = fu)

summary(aov(FAST_total_t2 ~ grupos, data = fu))

mean(bd$FAST_total_t2)
sd(bd$FAST_total_t2)
mean(mdd_eut$FAST_total_t2, na.rm = TRUE)
sd(mdd_eut$FAST_total_t2, na.rm = TRUE)
mean(mdd_cur$FAST_total_t2, na.rm = TRUE)
sd(mdd_cur$FAST_total_t2, na.rm = TRUE)

### COBRA ###

#quantile(bd$COBRA_soma_t2_16itens, probs = c(0.5,0.25,0.75))
#quantile(mdd_eut$COBRA_soma_t2_16itens, probs = c(0.5,0.25,0.75), na.rm = TRUE)
#quantile(mdd_cur$COBRA_soma_t2_16itens, probs = c(0.5,0.25,0.75), na.rm = TRUE)
#kruskal.test(COBRA_soma_t2_16itens ~ as.factor(grupos), data = fu)

summary(aov(COBRA_soma_t2_16itens ~ grupos, data = fu))

mean(bd$COBRA_soma_t2_16itens)
sd(bd$COBRA_soma_t2_16itens)
mean(mdd_eut$COBRA_soma_t2_16itens, na.rm = TRUE)
sd(mdd_eut$COBRA_soma_t2_16itens, na.rm = TRUE)
mean(mdd_cur$COBRA_soma_t2_16itens, na.rm = TRUE)
sd(mdd_cur$COBRA_soma_t2_16itens, na.rm = TRUE)

### LETTER-NUMBER SEQUENCING SUBTEST ###

#quantile(bd$NUMLETR_total, probs = c(0.5,0.25,0.75))
#quantile(mdd_eut$NUMLETR_total, probs = c(0.5,0.25,0.75), na.rm = TRUE)
#quantile(mdd_cur$NUMLETR_total, probs = c(0.5,0.25,0.75), na.rm = TRUE)
#kruskal.test(NUMLETR_total ~ as.factor(grupos), data = fu)

summary(aov(NUMLETR_total ~ grupos, data = fu))

mean(bd$NUMLETR_total)
sd(bd$NUMLETR_total)
mean(mdd_eut$NUMLETR_total)
sd(mdd_eut$NUMLETR_total)
mean(mdd_cur$NUMLETR_total)
sd(mdd_cur$NUMLETR_total)

#############################################################
###### TESTING FOR POSSIBLE CONFOUNDING VARIABLES IN BD GROUP
#############################################################

### YEARS OF EDUCATION ###
cor.test(bd$escolaridade, bd$COBRA_soma_t2_16itens, method="pearson")
cor.test(bd$escolaridade, bd$PSQI_total, method="pearson")
cor.test(bd$escolaridade, bd$FAST_total_t2, method="pearson")
cor.test(bd$escolaridade, bd$NUMLETR_total, method="pearson")

### CURRENT SUICIDE RISK ###
t.test(data = bd, COBRA_soma_t2_16itens ~ as.factor(riscodesuicidioatual))
t.test(data = bd, PSQI_total ~ as.factor(riscodesuicidioatual))
t.test(data = bd, FAST_total_t2 ~ as.factor(riscodesuicidioatual))
t.test(data = bd, NUMLETR_total ~ as.factor(riscodesuicidioatual))

### LIFETIME PSYCHIATRIC MEDICATION ###
t.test(data = bd, COBRA_soma_t2_16itens ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = bd, PSQI_total ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = bd, FAST_total_t2 ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = bd, NUMLETR_total ~ as.factor(lifetimepsychiatricmed_t2))

### HYPNOTICS/SEDATIVES ABUSE/DEPENDENCE ###
t.test(data = bd, COBRA_soma_t2_16itens ~ as.factor(hipnoticosabudep))
t.test(data = bd, PSQI_total ~ as.factor(hipnoticosabudep))
t.test(data = bd, FAST_total_t2 ~ as.factor(hipnoticosabudep))
t.test(data = bd, NUMLETR_total ~ as.factor(hipnoticosabudep))

### ANXIETY-RELATED COMORBID DISORDERS ###
t.test(data = bd, COBRA_soma_t2_16itens ~ as.factor(anxiety_dic))
t.test(data = bd, PSQI_total ~ as.factor(anxiety_dic))
t.test(data = bd, FAST_total_t2 ~ as.factor(anxiety_dic))
t.test(data = bd, NUMLETR_total ~ as.factor(anxiety_dic))

#######################################################################
###### TESTING FOR POSSIBLE CONFOUNDING VARIABLES IN EUTHYMIC MDD GROUP
#######################################################################

### YEARS OF EDUCATION ###
cor.test(mdd_eut$escolaridade, mdd_eut$COBRA_soma_t2_16itens, method="pearson")
cor.test(mdd_eut$escolaridade, mdd_eut$PSQI_total, method="pearson")
cor.test(mdd_eut$escolaridade, mdd_eut$FAST_total_t2, method="pearson")
cor.test(mdd_eut$escolaridade, mdd_eut$NUMLETR_total, method="pearson")

### CURRENT SUICIDE RISK ###
t.test(data = mdd_eut, COBRA_soma_t2_16itens ~ as.factor(riscodesuicidioatual))
t.test(data = mdd_eut, PSQI_total ~ as.factor(riscodesuicidioatual))
t.test(data = mdd_eut, FAST_total_t2 ~ as.factor(riscodesuicidioatual))
t.test(data = mdd_eut, NUMLETR_total ~ as.factor(riscodesuicidioatual))

### LIFETIME PSYCHIATRIC MEDICATION ###
t.test(data = mdd_eut, COBRA_soma_t2_16itens ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = mdd_eut, PSQI_total ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = mdd_eut, FAST_total_t2 ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = mdd_eut, NUMLETR_total ~ as.factor(lifetimepsychiatricmed_t2))

### HYPNOTICS/SEDATIVES ABUSE/DEPENDENCE ###
t.test(data = mdd_eut, COBRA_soma_t2_16itens ~ as.factor(hipnoticosabudep))
t.test(data = mdd_eut, PSQI_total ~ as.factor(hipnoticosabudep))
t.test(data = mdd_eut, FAST_total_t2 ~ as.factor(hipnoticosabudep))
t.test(data = mdd_eut, NUMLETR_total ~ as.factor(hipnoticosabudep))

### ANXIETY-RELATED COMORBID DISORDERS ###
t.test(data = mdd_eut, COBRA_soma_t2_16itens ~ as.factor(anxiety_dic))
t.test(data = mdd_eut, PSQI_total ~ as.factor(anxiety_dic))
t.test(data = mdd_eut, FAST_total_t2 ~ as.factor(anxiety_dic))
t.test(data = mdd_eut, NUMLETR_total ~ as.factor(anxiety_dic))

#################################################################################
###### TESTING FOR POSSIBLE CONFOUNDING VARIABLES IN MDD IN CURRENT EPISODE GROUP
#################################################################################

### YEARS OF EDUCATION ###
cor.test(mdd_cur$escolaridade, mdd_cur$COBRA_soma_t2_16itens, method="pearson")
cor.test(mdd_cur$escolaridade, mdd_cur$PSQI_total, method="pearson")
cor.test(mdd_cur$escolaridade, mdd_cur$FAST_total_t2, method="pearson")
cor.test(mdd_cur$escolaridade, mdd_cur$NUMLETR_total, method="pearson")

### CURRENT SUICIDE RISK ###
t.test(data = mdd_cur, COBRA_soma_t2_16itens ~ as.factor(riscodesuicidioatual))
t.test(data = mdd_cur, PSQI_total ~ as.factor(riscodesuicidioatual))
t.test(data = mdd_cur, FAST_total_t2 ~ as.factor(riscodesuicidioatual))
t.test(data = mdd_cur, NUMLETR_total ~ as.factor(riscodesuicidioatual))

### LIFETIME PSYCHIATRIC MEDICATION ###
t.test(data = mdd_cur, COBRA_soma_t2_16itens ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = mdd_cur, PSQI_total ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = mdd_cur, FAST_total_t2 ~ as.factor(lifetimepsychiatricmed_t2))
t.test(data = mdd_cur, NUMLETR_total ~ as.factor(lifetimepsychiatricmed_t2))

### HYPNOTICS/SEDATIVES ABUSE/DEPENDENCE ###
t.test(data = mdd_cur, COBRA_soma_t2_16itens ~ as.factor(hipnoticosabudep))
t.test(data = mdd_cur, PSQI_total ~ as.factor(hipnoticosabudep))
t.test(data = mdd_cur, FAST_total_t2 ~ as.factor(hipnoticosabudep))
t.test(data = mdd_cur, NUMLETR_total ~ as.factor(hipnoticosabudep))

### ANXIETY-RELATED COMORBID DISORDERS ###
t.test(data = mdd_cur, COBRA_soma_t2_16itens ~ as.factor(anxiety_dic))
t.test(data = mdd_cur, PSQI_total ~ as.factor(anxiety_dic))
t.test(data = mdd_cur, FAST_total_t2 ~ as.factor(anxiety_dic))
t.test(data = mdd_cur, NUMLETR_total ~ as.factor(anxiety_dic))

#################################################################################
###### TESTING FOR POSSIBLE CONFOUNDING VARIABLES WITH CORRELATION MATRIX
#################################################################################

####################################
###### CORRELATION TESTS IN BD GROUP
####################################

cor.test(bd$PSQI_total, bd$COBRA_soma_t2_16itens, method="pearson")
cor.test(bd$PSQI_total, bd$NUMLETR_total, method="pearson")
cor.test(bd$PSQI_total, bd$FAST_total_t2, method="pearson")

########################################################
###### CORRELATION TESTS IN MDD IN CURRENT EPISODE GROUP
########################################################

cor.test(mdd_cur$PSQI_total, mdd_cur$COBRA_soma_t2_16itens, method="pearson")
cor.test(mdd_cur$PSQI_total, mdd_cur$NUMLETR_total, method="pearson")
cor.test(mdd_cur$PSQI_total, mdd_cur$FAST_total_t2, method="pearson")

########################################
###### CORRELATION TESTS IN EUTHYMIC MDD
########################################

cor.test(mdd_eut$PSQI_total, mdd_eut$COBRA_soma_t2_16itens, method="pearson")
cor.test(mdd_eut$PSQI_total, mdd_eut$NUMLETR_total, method="pearson")
cor.test(mdd_eut$PSQI_total, mdd_eut$FAST_total_t2, method="pearson")

##################################################
###### KRUSKAL-WALLIS AND POST-HOCS BETWEEN GROUPS
##################################################

# 0 = BD; 1 = MDD in current episode; 2 = Euthymic MDD

#kruskal.test(PSQI_total ~ as.factor(grupos), data = fu)
#dunnTest(PSQI_total ~ as.factor(grupos), data = fu, method = "holm")
#dscfAllPairsTest(PSQI_total ~ as.factor(grupos), data=fu)
#
#kruskal.test(FAST_total_t2 ~ as.factor(grupos), data = fu)
#dunnTest(FAST_total_t2 ~ as.factor(grupos), data = fu, method = "holm")
#dscfAllPairsTest(FAST_total_t2 ~ as.factor(grupos), data = fu)
#
#kruskal.test(NUMLETR_total ~ as.factor(grupos), data = fu)
#dunnTest(NUMLETR_total ~ as.factor(grupos), data = fu, method = "holm")
#dscfAllPairsTest(NUMLETR_total ~ as.factor(grupos), data = fu)
#
#kruskal.test(COBRA_soma_t2_16itens ~ as.factor(grupos), data = fu)
#dunnTest(COBRA_soma_t2_16itens ~ as.factor(grupos), data = fu, method = "holm")
#dscfAllPairsTest(COBRA_soma_t2_16itens ~ as.factor(grupos), data = fu)

##################################################
###### ANOVA AND TUKEY POST-HOCS BETWEEN GROUPS
##################################################

summary(aov(PSQI_total ~ grupos, data = fu))
TukeyHSD(aov(PSQI_total ~ grupos, data = fu))

summary(aov(FAST_total_t2 ~ grupos, data = fu))
TukeyHSD(aov(FAST_total_t2 ~ grupos, data = fu))

summary(aov(COBRA_soma_t2_16itens ~ grupos, data = fu))
TukeyHSD(aov(COBRA_soma_t2_16itens ~ grupos, data = fu))

summary(aov(NUMLETR_total ~ grupos, data = fu))
TukeyHSD(aov(NUMLETR_total ~ grupos, data = fu))

###############################
###### CRUDE LINEAR REGRESSIONS
###############################

# BD group
summary(fast_bd <- lm(FAST_total_t2 ~ PSQI_total, data = bd))
confint(fast_bd)
summary(cobra_bd <- lm(COBRA_soma_t2_16itens ~ PSQI_total, data = bd))
confint(cobra_bd)
summary(numletr_bd <- lm(NUMLETR_total ~ PSQI_total, data = bd))
confint(numletr_bd)

# Euthymic MDD group
summary(fast_mdd_eut <- lm(FAST_total_t2 ~ PSQI_total, data = mdd_eut))
confint(fast_mdd_eut)
summary(cobra_mdd_eut <-lm(COBRA_soma_t2_16itens ~ PSQI_total, data = mdd_eut))
confint(cobra_mdd_eut)
summary(numletr_mdd_eut <- lm(NUMLETR_total ~ PSQI_total, data = mdd_eut))
confint(numletr_mdd_eut)

# MDD in current episode group
summary(fast_mdd_cur <- lm(FAST_total_t2 ~ PSQI_total, data = mdd_cur))
confint(fast_mdd_cur)
summary(cobra_mdd_cur <- lm(COBRA_soma_t2_16itens ~ PSQI_total, data = mdd_cur))
confint(cobra_mdd_cur)
summary(numletr_mdd_cur <- lm(NUMLETR_total ~ PSQI_total, data = mdd_cur))
confint(numletr_mdd_cur)

##############################################################
###### LINEAR REGRESSIONS ACCOUNTING FOR CONFOUNDING VARIABLES
##############################################################

# BD group
summary(fast_bd_adj <- lm(FAST_total_t2 ~ PSQI_total + escolaridade + lifetimepsychiatricmed_t2 + hipnoticosabudep, data = bd))
confint(fast_bd_adj)
summary(cobra_bd_adj <- lm(COBRA_soma_t2_16itens ~ PSQI_total + escolaridade + lifetimepsychiatricmed_t2, data = bd))
confint(cobra_bd_adj)
summary(numletr_bd_adj <- lm(NUMLETR_total ~ PSQI_total + escolaridade + hipnoticosabudep, data = bd))
confint(numletr_bd_adj)

# Euthymic MDD group
summary(fast_mdd_eut_adj <- lm(FAST_total_t2 ~ PSQI_total + anxiety_dic, data = mdd_eut))
confint(fast_mdd_eut_adj)
summary(cobra_mdd_eut_adj <-lm(COBRA_soma_t2_16itens ~ PSQI_total + anxiety_dic, data = mdd_eut))
confint(cobra_mdd_eut_adj)
summary(numletr_mdd_eut_adj <- lm(NUMLETR_total ~ PSQI_total, data = mdd_eut))
confint(numletr_mdd_eut_adj)

# MDD in current episode group
summary(fast_mdd_cur_adj <- lm(FAST_total_t2 ~ PSQI_total + escolaridade + riscodesuicidioatual + anxiety_dic, data = mdd_cur))
confint(fast_mdd_cur_adj)
summary(cobra_mdd_cur_adj <- lm(COBRA_soma_t2_16itens ~ PSQI_total + escolaridade + riscodesuicidioatual, data = mdd_cur))
confint(cobra_mdd_cur_adj)
summary(numletr_mdd_cur_adj <- lm(NUMLETR_total ~ PSQI_total + escolaridade + riscodesuicidioatual, data = mdd_cur))
confint(numletr_mdd_cur_adj)
