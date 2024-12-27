import numpy as np
import pandas as pd

df = pd.read_csv('dataset_preprocessed.csv')


rows = [
'N',
'Age, mean (SD), y',
'Race & Ethnicity, n(%)',
'    Asian',
'    Black',
'    Hispanic',
'    White',
'    Other',
'Education, n(%)',
'    Some elementary school',
'    Elementary school',
'    Some high school',
'    High school',
'    Some college',
'    College',
'    Some graduate school',
'    Graduate school',
'Body mass index, median (IQR), kg/m2',
'Apnea-hypopnea index, median (IQR), /hour',
'APOE e4 carrier, n(%)',
'Diagnosis, n(%)',
'    Hypertension',
'    Type 2 diabetes',
'    Stroke',
'Sleep medication use, n(%)',
'    Benzodiazepine',
'    Antidepressant',
'    Zolpidem',
'    Opiod',
'    Anti-inflammation',
'Teng Modified Mini-Mental State (3MS) Test score at Visit 2, median (IQR)',
]


col = 'Value'
df_res = pd.DataFrame(data={'Name':rows})
df_res[col] = np.nan

df_res.loc[rows.index('N'), col] = len(df)
df_res.loc[rows.index('Age, mean (SD), y'), col] = f'{df.Age_VS1.mean():.1f} ({df.Age_VS1.std():.1f})'

df_res.loc[rows.index('    Asian'), col] = f'{df.Race_Asian.sum()} ({df.Race_Asian.mean()*100:.1f}%)'
df_res.loc[rows.index('    Black'), col] = f'{df.Race_Black.sum()} ({df.Race_Black.mean()*100:.1f}%)'
df_res.loc[rows.index('    Hispanic'), col] = f'{df.Race_Hispanic.sum()} ({df.Race_Hispanic.mean()*100:.1f}%)'
n = len(df)-df.Race_Asian.sum()-df.Race_Black.sum()-df.Race_Hispanic.sum()-df.Race_Other.sum()
df_res.loc[rows.index('    White'), col] = f'{n} ({n/len(df)*100:.1f}%)'
df_res.loc[rows.index('    Other'), col] = f'{df.Race_Other.sum()} ({df.Race_Other.mean()*100:.1f}%)'

df_res.loc[rows.index('    Some elementary school'), col] = f'{(df.Educ==1).sum()} ({(df.Educ==1).mean()*100:.1f}%)'
df_res.loc[rows.index('    Elementary school'), col] = f'{(df.Educ==2).sum()} ({(df.Educ==2).mean()*100:.1f}%)'
df_res.loc[rows.index('    Some high school'), col] = f'{(df.Educ==3).sum()} ({(df.Educ==3).mean()*100:.1f}%)'
df_res.loc[rows.index('    High school'), col] = f'{(df.Educ==4).sum()} ({(df.Educ==4).mean()*100:.1f}%)'
df_res.loc[rows.index('    Some college'), col] = f'{(df.Educ==5).sum()} ({(df.Educ==5).mean()*100:.1f}%)'
df_res.loc[rows.index('    College'), col] = f'{(df.Educ==6).sum()} ({(df.Educ==6).mean()*100:.1f}%)'
df_res.loc[rows.index('    Some graduate school'), col] = f'{(df.Educ==7).sum()} ({(df.Educ==7).mean()*100:.1f}%)'
df_res.loc[rows.index('    Graduate school'), col] = f'{(df.Educ==8).sum()} ({(df.Educ==8).mean()*100:.1f}%)'


q1, q2, q3 = np.nanpercentile(df.BMI_VS1, (25,50,75))
df_res.loc[rows.index('Body mass index, median (IQR), kg/m2'), col] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.AHI4, (25,50,75))
df_res.loc[rows.index('Apnea-hypopnea index, median (IQR), /hour'), col] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.loc[rows.index('APOE e4 carrier, n(%)'), col] = f'{(df.APOE4Count>0).sum():.0f} ({(df.APOE4Count>0).mean()*100:.1f}%)'

df_res.loc[rows.index('    Hypertension'), col] = f'{df.MH_HTN_VS1.sum():.0f} ({df.MH_HTN_VS1.mean()*100:.1f}%)'
df_res.loc[rows.index('    Type 2 diabetes'), col] = f'{df.MH_DB2_VS1.sum():.0f} ({df.MH_DB2_VS1.mean()*100:.1f}%)'
df_res.loc[rows.index('    Stroke'), col] = f'{df.MH_Stroke_VS1.sum():.0f} ({df.MH_Stroke_VS1.mean()*100:.1f}%)'

df_res.loc[rows.index('    Benzodiazepine'), col] = f'{df.Med_Benzo_VS1.sum():.0f} ({df.Med_Benzo_VS1.mean()*100:.1f}%)'
df_res.loc[rows.index('    Antidepressant'), col] = f'{df.Med_Antidep_VS1.sum():.0f} ({df.Med_Antidep_VS1.mean()*100:.1f}%)'
df_res.loc[rows.index('    Zolpidem'), col] = f'{df.Med_Zolpidem_VS1.sum():.0f} ({df.Med_Zolpidem_VS1.mean()*100:.1f}%)'
df_res.loc[rows.index('    Opiod'), col] = f'{df.Med_Opiod_VS1.sum():.0f} ({df.Med_Opiod_VS1.mean()*100:.1f}%)'
df_res.loc[rows.index('    Anti-inflammation'), col] = f'{df.Med_AntiInfl_VS1.sum():.0f} ({df.Med_AntiInfl_VS1.mean()*100:.1f}%)'

q1, q2, q3 = np.nanpercentile(df.Teng3MSScore_V2, (25,50,75))
df_res.loc[rows.index('Teng Modified Mini-Mental State (3MS) Test score at Visit 2, median (IQR)'), col] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'

print(df_res)
df_res.to_excel('table1.xlsx', index=False)

