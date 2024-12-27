library(readxl)
library(mediation)

# load dataset
df.all <- read.csv('/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/dataset_preprocessed.csv')

covariates.basic <- c('Educ', 'Race_Black','Race_Asian', 'Race_Hispanic', 'Race_Other', "APOE4Count")
#covariates.v1 <- c('Age_V1', 'BMI_V1', 'Smoking', 'Alcohol', 'Med_AntiInfl_V1', "MH_HTN_V1", "MH_Stroke_V1", "MH_DB2_V1")
covariates.vs1 <- c("Age_VS1", "BMI_VS1", "MH_HTN_VS1", "MH_Stroke_VS1", "MH_DB2_VS1", "Med_Benzo_VS1", "Med_Antidep_VS1", "Med_Zolpidem_VS1", "Med_Opiod_VS1", "Med_AntiInfl_VS1", 'AHI4')
covariates.basic.str <- paste0(covariates.basic, collapse='+')
covariates.vs1.str <- paste0(covariates.vs1, collapse='+')

med.res.path <- '/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/mediation_results_exposure_biomarkervs1_mediator_sleepvs1_outcomev2_refined.xlsx'
df.med <- read_excel(med.res.path)
df.med <- df.med[df.med[,'Significant(ie)']==TRUE,]

sims <- 10000
ii <- 0
sens.res.all <- list()

for (i in 1:nrow(df.med)) {
  exposure <- df.med$exposure[i]
  mediator <- df.med$mediator[i]
  outcome <- df.med$outcome[i]
  print(sprintf('[%d/%d] exposure = %s, mediator = %s, outcome = %s', i, nrow(df.med), exposure, mediator, outcome))
  
  df <- data.frame(df.all)
  #df[,exposure] <- cut(df[,exposure], quantile(df[,exposure], c(0,0.5,1)), include.lowest=TRUE, labels=FALSE, right=TRUE)-1
  
  # outcome model: P(outcome | exposure, mediator, covariates)
  #TODO exposure-mediator interaction, but exposure is converted to 1,2,3, not
  formula.y.str <- sprintf('%s ~ %s + %s + %s + %s', outcome, exposure, mediator, covariates.basic.str, covariates.vs1.str)#, covariates.v1.str
  formula.y <- as.formula(formula.y.str)
  model.y <- lm(formula.y, df)
  
  # mediator model: P(mediator | exposure, covariates)
  formula.m.str <- sprintf('%s ~ %s + %s + %s', mediator, exposure, covariates.basic.str, covariates.vs1.str)
  formula.m <- as.formula(formula.m.str)
  model.m <- lm(formula.m, df)
  
  low <- quantile(df[,exposure], 0.25, na.rm=TRUE)
  high <- quantile(df[,exposure], 0.75, na.rm=TRUE)
  #low <- 0
  #high <- 1
  set.seed(2024)
  med.res <- mediate(model.m, model.y,
                     treat=exposure, mediator=mediator,
                     covariates = c(covariates.basic, covariates.vs1),
                     control.value = low, treat.value = high,
                     boot=TRUE, sims=sims, boot.ci.type='bca')
  #res.summary <- summary(med.res)
  #print(res.summary)
  
  #med.int <- test.TMint(med.res)
  sens.res <- medsens(med.res, rho.by=0.01, effect.type='indirect', sims=1000)
  sens.res.all[[i]] <-sens.res
  
  med.mod.age <- test.modmed(med.res,covariates.1=list(Age_VS1=70),covariates.2=list(Age_VS1=80),sims=sims)
  sink(sprintf("/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/sens-ie-age-%s-%s-%s.txt", exposure, mediator, outcome))
  print(med.mod.age)
  sink()
  
  med.mod.bmi <- test.modmed(med.res,covariates.1=list(BMI_VS1=25),covariates.2=list(BMI_VS1=30),sims=sims)
  sink(sprintf("/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/sens-ie-BMI-%s-%s-%s.txt", exposure, mediator, outcome))
  print(med.mod.bmi)
  sink()
}
saveRDS(sens.res.all, "/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/sens.res.all_exposure_biomarkervs1_mediator_sleepvs1_outcomev2.rds")
