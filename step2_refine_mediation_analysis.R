library(readxl)
library(mediation)

# load dataset
df.all <- read.csv('/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/dataset_preprocessed.csv')

covariates.basic <- c('Educ', 'Race_Black','Race_Asian', 'Race_Hispanic', 'Race_Other', "APOE4Count")
#covariates.v1 <- c('Age_V1', 'BMI_V1', 'Smoking', 'Alcohol', 'Med_AntiInfl_V1', "MH_HTN_V1", "MH_Stroke_V1", "MH_DB2_V1")
covariates.vs1 <- c("Age_VS1", "BMI_VS1", "MH_HTN_VS1", "MH_Stroke_VS1", "MH_DB2_VS1", "Med_Benzo_VS1", "Med_Antidep_VS1", "Med_Zolpidem_VS1", "Med_Opiod_VS1", "Med_AntiInfl_VS1", 'AHI4')
covariates.basic.str <- paste0(covariates.basic, collapse='+')
covariates.vs1.str <- paste0(covariates.vs1, collapse='+')


med.res.path <- '/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/mediation_results_exposure_biomarkervs1_mediator_sleepvs1_outcomev2.csv'
#med.res.path <- '/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/mediation_results_exposure_sleepvs1_mediator_biomarkervs1_outcomev2.csv'
df.med <- read.csv(med.res.path)
df.med <- df.med[order(df.med$ie.avg.p),]
df.med <- df.med[(df.med$ie.avg.p<0.1)|(df.med$de.avg.p<0.05),]

sims <- 10000
df.res <- list()
ii <- 0
result.path <- 'mediation_results_exposure_biomarkervs1_mediator_sleepvs1_outcomev2_refined2.csv'
#result.path <- 'mediation_results_exposure_sleepvs1_mediator_biomarkervs1_outcomev2_refined2.csv'

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
  res.summary <- summary(med.res)
  print(res.summary)
  
  model.m.coef <- coefficients(res.summary$model.m)
  model.m.coef <- model.m.coef[2:length(model.m.coef)]
  names(model.m.coef)[names(model.m.coef)==exposure] <- 'exposure'
  
  model.m.ci <- confint(res.summary$model.m)
  model.m.ci <- model.m.ci[2:nrow(model.m.ci),]
  rownames(model.m.ci)[rownames(model.m.ci)==exposure] <- 'exposure'
  
  model.m.p <- coefficients(summary(res.summary$model.m))[,'Pr(>|t|)']
  model.m.p <- model.m.p[2:length(model.m.p)]
  names(model.m.p)[names(model.m.p)==exposure] <- 'exposure'
  
  model.y.coef <- coefficients(res.summary$model.y)
  model.y.coef <- model.y.coef[2:length(model.y.coef)]
  names(model.y.coef)[names(model.y.coef)==exposure] <- 'exposure'
  names(model.y.coef)[names(model.y.coef)==mediator] <- 'mediator'
  
  model.y.ci <- confint(res.summary$model.y)
  model.y.ci <- model.y.ci[2:nrow(model.y.ci),]
  rownames(model.y.ci)[rownames(model.y.ci)==exposure] <- 'exposure'
  rownames(model.y.ci)[rownames(model.y.ci)==mediator] <- 'mediator'
  
  model.y.p <- coefficients(summary(res.summary$model.y))[,'Pr(>|t|)']
  model.y.p <- model.y.p[2:length(model.y.p)]
  names(model.y.p)[names(model.y.p)==exposure] <- 'exposure'
  names(model.y.p)[names(model.y.p)==mediator] <- 'mediator'
  
  names(res.summary$d.avg.ci) <- c('2.5%', '97.5%')
  names(res.summary$z.avg.ci) <- c('2.5%', '97.5%')
  names(res.summary$n.avg.ci) <- c('2.5%', '97.5%')
  names(res.summary$tau.ci) <- c('2.5%', '97.5%')
  names(res.summary$d0.ci) <- c('2.5%', '97.5%')
  names(res.summary$d1.ci) <- c('2.5%', '97.5%')
  names(res.summary$z0.ci) <- c('2.5%', '97.5%')
  names(res.summary$z1.ci) <- c('2.5%', '97.5%')
  names(res.summary$n0.ci) <- c('2.5%', '97.5%')
  names(res.summary$n1.ci) <- c('2.5%', '97.5%')
  
  res <- c(nobs=res.summary$nobs,
           exposure0=as.numeric(res.summary$control.value), exposure1=as.numeric(res.summary$treat.value),
           ie.avg=res.summary$d.avg, ie.avg.ci=res.summary$d.avg.ci, ie.avg.p=res.summary$d.avg.p,
           de.avg=res.summary$z.avg, de.avg.ci=res.summary$z.avg.ci, de.avg.p=res.summary$z.avg.p,
           ip.avg=res.summary$n.avg, ip.avg.ci=res.summary$n.avg.ci, ip.avg.p=res.summary$n.avg.p,
           te=res.summary$tau.coef, te.ci=res.summary$tau.ci, te.p=res.summary$tau.p,
           ie0=res.summary$d0, ie0.ci=res.summary$d0.ci, ie0.p=res.summary$d0.p,
           ie1=res.summary$d1, ie1.ci=res.summary$d1.ci, ie1.p=res.summary$d1.p,
           de0=res.summary$z0, de0.ci=res.summary$z0.ci, de0.p=res.summary$z0.p,
           de1=res.summary$z1, de1.ci=res.summary$z1.ci, de1.p=res.summary$z1.p,
           ip0=res.summary$n0, ip0.ci=res.summary$n0.ci, ip0.p=res.summary$n0.p,
           ip1=res.summary$n1, ip1.ci=res.summary$n1.ci, ip1.p=res.summary$n1.p,
           model.m.coef=model.m.coef, model.m.ci.lb=model.m.ci[,1], model.m.ci.ub=model.m.ci[,2], model.m.p=model.m.p,
           model.y.coef=model.y.coef, model.y.ci.lb=model.y.ci[,1], model.y.ci.ub=model.y.ci[,2], model.y.p=model.y.p
  )
  #boot=res.summary$boot, boot.ci.type=res.summary$boot.ci.type, INT=res.summary$INT, nsim=res.summary$sims
  res <- as.data.frame(t(res))
  res <- cbind(exposure=res.summary$treat, mediator=res.summary$mediator, outcome=outcome, res)
  ii <- ii+1
  df.res[[ii]] <- res
  
  df.res2 <- do.call(rbind, df.res)
  write.csv(df.res2, file=result.path, row.names=F, quote=F, na="")
  
  #med.int <- test.TMint(med.res)
  #med.mod.age <- test.modmed(med.res,covariates.1=list(Age_VS1=70),covariates.2=list(Age_VS1=90),sims=sims)
  #med.mod.bmi <- test.modmed(med.res,covariates.1=list(BMI_VS1=25),covariates.2=list(BMI_VS1=40),sims=sims)
  #sens.res <- medsens(med.res, rho.by=0.01, effect.type='both', sims=10)
  #plot(sens.res, sens.par='rho')
  #plot(sens.res, sens.par='R2')
}

df.res2 <- do.call(rbind, df.res)
print(df.res2)
write.csv(df.res2, file=result.path, row.names=F, quote=F, na="")
