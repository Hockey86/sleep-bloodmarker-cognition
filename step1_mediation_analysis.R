library(readxl)
library(mediation)

# load dataset
data.path <- '/data/haoqisun/inflammation-sleep-dementia/dataset_MrOS.xlsx'
#data1 <- read_excel(data.path, sheet='Biomarkers_V1')
data2 <- read_excel(data.path, sheet='Biomarkers_VS1')
data3 <- read_excel(data.path, sheet='Covariates_V1')
data4 <- read_excel(data.path, sheet='Covariates_VS1')
data5 <- read_excel(data.path, sheet='Covariates_V2')
data6 <- read_excel(data.path, sheet='Sleep_VS1')

data <- cbind(data2, data3[,2:ncol(data3)], data4[,2:ncol(data4)], data5[,2:ncol(data5)], data6[,2:ncol(data6)])

#data$Race_Other2 <- as.integer((data$Race_Hispanic+data$Race_Other)>0)

# log-transformation to the cytokine levels
cytokine.names <- c(names(data2)[2:ncol(data2)])
data[,cytokine.names] <- log(data[,cytokine.names])

eeg.names <- c('delta_dbs_N3_C', 'delta_rel_N3_C', 'alpha_dbs_N3_C', 'alpha_rel_N3_C', 'delta_dbs_N2_C', 'delta_rel_N2_C', 'theta_dbs_N1_C', 'theta_rel_N1_C', 'delta_dbs_R_C', 'delta_rel_R_C', 'theta_dbs_R_C', 'theta_rel_R_C', 'delta_slope_N2N3_C', 'SP_FFT_all_C', 'SP_AMP_all_C', 'SP_CDENS_all_C', 'SP_CHIRP_all_C', 'SP_COUPL_MAG_all_C', 'SP_COUPL_OVERLAP_all_C', 'SP_DENS_all_C', 'SP_ISA_S_all_C', 'SP_AMP_slow_C', 'SP_CDENS_slow_C', 'SP_CHIRP_slow_C', 'SP_COUPL_MAG_slow_C', 'SP_COUPL_OVERLAP_slow_C', 'SP_DENS_slow_C', 'SP_ISA_S_slow_C', 'SP_AMP_fast_C', 'SP_CDENS_fast_C', 'SP_CHIRP_fast_C', 'SP_COUPL_MAG_fast_C', 'SP_COUPL_OVERLAP_fast_C', 'SP_DENS_fast_C', 'SP_ISA_S_fast_C', 'SO_SLOPE_NEG1_C', 'SO_SLOPE_POS1_C', 'SO_RATE_C', 'SO_SLOPE_C', 'SO_DUR_C', 'SO_P2P_C')

#result.path <- 'mediation_results_exposure_biomarkervs1_mediator_sleepvs1_outcomev2.csv'
result.path <- 'mediation_results_exposure_sleepvs1_mediator_biomarkervs1_outcomev2.csv'

exposures <- eeg.names#[1:2]
mediators <- cytokine.names#[1:2]
outcomes <- c("Teng3MSScore_V2")#, "TrialsBTime_V2")

# define outcome model and mediator model
#    outcome: dementia
#    exposure: cytokine
#    mediator: sleep pattern

covariates.basic <- c('Educ', 'Race_Black','Race_Asian', 'Race_Hispanic', 'Race_Other', "APOE4Count")
#covariates.v1 <- c('Age_V1', 'BMI_V1', 'Smoking', 'Alcohol', 'Med_AntiInfl_V1', "MH_HTN_V1", "MH_Stroke_V1", "MH_DB2_V1")
covariates.vs1 <- c("Age_VS1", "BMI_VS1", "MH_HTN_VS1", "MH_Stroke_VS1", "MH_DB2_VS1", "Med_Benzo_VS1", "Med_Antidep_VS1", "Med_Zolpidem_VS1", "Med_Opiod_VS1", "Med_AntiInfl_VS1", 'AHI4','AHI3')

# deal with missing value
# we can use the simplest approach: exclude any subjects with any missing value
cols <- c('ID', exposures,  mediators,outcomes,covariates.basic, covariates.vs1)#, covariates.v1
print(dim(data[,cols]))
ids <- complete.cases(data[,cols])
data <- data[ids,cols]
print(dim(data))

write.csv(data, file='/data/haoqisun/inflammation-sleep-dementia/mediation-analysis/dataset_preprocessed.csv', row.names=F, quote=F)

covariates.basic.str <- paste0(covariates.basic, collapse='+')
#covariates.v1.str <- paste0(covariates.v1, collapse='+')
covariates.vs1.str <- paste0(covariates.vs1, collapse='+')

df.res <- list()
ii <- 0
for (ei in 1:length(exposures)) {
  exposure <- exposures[ei]
  for (mi in 1:length(mediators)) {
    mediator <- mediators[mi]
    for (oi in 1:length(outcomes)) {
      outcome <- outcomes[oi]
      print(sprintf('[%d/%d] exposure = %s, [%d/%d] mediator = %s, [%d/%d] outcome = %s', ei, length(exposures), exposure, mi, length(mediators), mediator, oi, length(outcomes), outcome))
      #data <- data.frame(data.all)
      
      # outcome model: P(outcome | exposure, mediator, covariates)
      formula.y.str <- sprintf('%s ~ %s + %s + %s + %s', outcome, exposure, mediator, covariates.basic.str, covariates.vs1.str)#, covariates.v1.str
      formula.y <- as.formula(formula.y.str)
      model.y <- lm(formula.y, data)
      
      # mediator model: P(mediator | exposure, covariates)
      formula.m.str <- sprintf('%s ~ %s + %s + %s', mediator, exposure, covariates.basic.str, covariates.vs1.str)
      formula.m <- as.formula(formula.m.str)
      model.m <- lm(formula.m, data)
      
      # define the low and high levels for the exposure, since our exposure is not binary, it is a continuous number
      
      # to start simple, we can set the low level as quartile 1 (Q1)
      low <- quantile(data[,exposure], 0.25, na.rm=TRUE)
      high <- quantile(data[,exposure], 0.75, na.rm=TRUE)
      
      med.res <- mediate(model.m, model.y,
                        treat=exposure, mediator=mediator,
                        covariates = c(covariates.basic, covariates.vs1), outcome = outcome,
                        control.value = low, treat.value = high,
                        robustSE=TRUE)#boot=TRUE, 
      
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
    } 
  }
}

df.res2 <- do.call(rbind, df.res)
print(df.res2)
write.csv(df.res2, file=result.path, row.names=F, quote=F, na="")

#TODO nonlinear model?
#TODO different levels
