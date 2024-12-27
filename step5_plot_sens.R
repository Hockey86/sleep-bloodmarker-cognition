library(readxl)
library(mediation)

data.dir <- '/data/haoqisun/inflammation-sleep-dementia/mediation-analysis'


# IE for exposure=biomarkervs1

res <- readRDS(file.path(data.dir, "sens-analysis/exposure_biomarkervs1/sens.res.all_exposure_biomarkervs1_mediator_sleepvs1_outcomev2.rds"))
med.res.path <- file.path(data.dir, 'mediation_results_exposure_biomarkervs1_mediator_sleepvs1_outcomev2_refined.xlsx')
df.med <- read_excel(med.res.path)
df.med <- df.med[df.med[,'Significant(ie)']==T,]

stopifnot(nrow(df.med)==length(res))

for (i in 1:nrow(df.med)) {
  exposure <- df.med$exposure[i]
  mediator <- df.med$mediator[i]
  outcome <- df.med$outcome[i]
  sens.res <- res[[i]]
  
  sink(file.path(data.dir, sprintf("sens-analysis/exposure_biomarkervs1/sens-ie-%s-%s-%s.txt", exposure, mediator, outcome)))
  print(summary(sens.res))
  sink()
  
  #lower <- quantile((sens.res$lower.d1+sens.res$lower.d0)/2,0.05)
  #upper <- quantile((sens.res$upper.d1+sens.res$upper.d0)/2,0.95)
  ylim <- 0.5#max(abs(lower), abs(upper))
  ylim <- c(-ylim, ylim)
  
  png(file.path(data.dir, sprintf("sens-analysis/exposure_biomarkervs1/sens-ie-R2-%s-%s-%s.png", exposure, mediator, outcome)), units="in", width=4, height=4, res=300)
  par(mar=c(4,4,1,1))
  plot(sens.res, sens.par='R2', main=NA)
  dev.off()
  
  png(file.path(data.dir, sprintf("sens-analysis/exposure_biomarkervs1/sens-ie-rho-%s-%s-%s.png", exposure, mediator, outcome)), units="in", width=5, height=5/3*2, res=300)
  par(mar=c(4,4,1,1))
  plot(sens.res, sens.par='rho', ylim=ylim, main=NA)
  dev.off()
}



# DE for exposure=sleepvs1
res <- readRDS(file.path(data.dir, "sens-analysis/exposure_sleepvs1/sens.res.all_exposure_sleepvs1_mediator_biomarkervs1_outcomev2.rds"))
med.res.path <- file.path(data.dir, 'mediation_results_exposure_sleepvs1_mediator_biomarkervs1_outcomev2_refined.xlsx')
df.med <- read_excel(med.res.path, sheet='Significant_direct_effect')

stopifnot(nrow(df.med)==length(res))

for (i in 1:nrow(df.med)) {
  exposure <- df.med$exposure[i]
  outcome <- df.med$outcome[i]
  sens.res <- res[[i]]
  
  sink(file.path(data.dir, sprintf("sens-analysis/exposure_sleepvs1/sens-de-%s-%s.txt", exposure, outcome)))
  print(summary(sens.res))
  sink()
  
  lower <- quantile((sens.res$lower.z1+sens.res$lower.z0)/2,0.05)
  upper <- quantile((sens.res$upper.z1+sens.res$upper.z0)/2,0.95)
  ylim <- max(abs(lower), abs(upper))
  ylim <- c(-ylim, ylim)
  
  png(file.path(data.dir, sprintf("sens-analysis/exposure_sleepvs1/sens-de-R2-%s-%s.png", exposure, outcome)), units="in", width=4, height=4, res=300)
  par(mar=c(4,4,1,1))
  plot(sens.res, sens.par='R2', main=NA)
  dev.off()
  
  png(file.path(data.dir, sprintf("sens-analysis/exposure_sleepvs1/sens-de-rho-%s-%s.png", exposure, outcome)), units="in", width=5, height=5/3*2, res=300)
  par(mar=c(4,4,1,1))
  plot(sens.res, sens.par='rho', ylim=ylim, main=NA)
  dev.off()
}