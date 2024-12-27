import pandas as pd
import pingouin


df=pd.read_csv('dataset_preprocessed.csv')

covars = ['Educ', 'Race_Black','Race_Asian', 'Race_Hispanic', 'Race_Other', "APOE4Count","Age_VS1", "BMI_VS1", "MH_HTN_VS1", "MH_Stroke_VS1", "MH_DB2_VS1", "Med_Benzo_VS1", "Med_Antidep_VS1", "Med_Zolpidem_VS1", "Med_Opiod_VS1", "Med_AntiInfl_VS1", 'AHI4']

variables = [
['Leptin_VS1', 'SP_DENS_fast_C'],
['SP_DENS_fast_C', 'Teng3MSScore_V2'],
['Leptin_VS1', 'Teng3MSScore_V2'],
['Leptin_VS1', 'SP_COUPL_OVERLAP_all_C'],
['SP_COUPL_OVERLAP_all_C', 'Teng3MSScore_V2'],
['theta_dbs_N1_C', 'Teng3MSScore_V2'],
['theta_rel_N1_C', 'Teng3MSScore_V2'],
['theta_dbs_R_C', 'Teng3MSScore_V2'],
['SP_DENS_all_C', 'Teng3MSScore_V2'],
['SP_DENS_fast_C', 'Teng3MSScore_V2'],
['SP_DENS_slow_C', 'Teng3MSScore_V2'],
['SP_CDENS_all_C', 'Teng3MSScore_V2'],
['SP_CDENS_fast_C', 'Teng3MSScore_V2'],
['SP_COUPL_OVERLAP_all_C', 'Teng3MSScore_V2'],
['SP_COUPL_OVERLAP_slow_C', 'Teng3MSScore_V2'],
#['SP_COUPL_MAG_fast_C', 'Teng3MSScore_V2'],
]

df_res = []
for v in variables:
    df_ = pingouin.partial_corr(df,x=v[0], y=v[1], covar=covars, method='spearman')
    df_.insert(0, 'Variable2', v[1])
    df_.insert(0, 'Variable1', v[0])
    df_res.append(df_)
df_res = pd.concat(df_res, axis=0, ignore_index=True)

print(df_res)
df_res.to_excel('corrs_exposure_mediator_outcome.xlsx', index=False)


