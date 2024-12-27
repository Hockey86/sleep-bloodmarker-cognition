import numpy as np
import pandas as pd


def get_bmi_cat(bmi):
    if bmi<18.5:
        return 0
    elif bmi<35:
        return 1
    else:
        return 2
        

def get_ahi_cat(ahi):
    if ahi<5:
        return 0
    elif ahi<15:
        return 1
    elif ahi<30:
        return 2
    else:
        return 3
        
        
def main():
    # load dataset
    df_path =  '/data/haoqisun/inflammation-sleep-dementia/dataset_MrOS.xlsx'
    with pd.ExcelFile(df_path) as f:
        #df1 = pd.read_excel(f, 'Biomarkers_V1')
        df2 = pd.read_excel(f, 'Biomarkers_VS1')
        #df3 = pd.read_excel(f, 'Covariates_V1')
        df4 = pd.read_excel(f, 'Covariates_VS1')
        df5 = pd.read_excel(f, 'Covariates_V2')
        df6 = pd.read_excel(f, 'PSG_VS1')
        df7 = pd.read_excel(f, 'SleepEEG_VS1')

    df = pd.concat([df1, df2.iloc[:,1:], df3.iloc[:,1:], df4.iloc[:,1:], df5.iloc[:,1:], df6.iloc[:,1:], df7.iloc[:,1:]], axis=1)

    # log-transformation to the cytokine levels
    cytokine_names = list(df1.columns[1:])+list(df2.columns[1:])
    df.loc[:,cytokine_names] = np.log(df.loc[:,cytokine_names])

    eeg_names = ['delta_dbs_N3_C', 'delta_rel_N3_C', 'alpha_dbs_N3_C', 'alpha_rel_N3_C', 'delta_dbs_N2_C', 'delta_rel_N2_C', 'theta_dbs_N1_C', 'theta_rel_N1_C', 'delta_dbs_R_C', 'delta_rel_R_C', 'theta_dbs_R_C', 'theta_rel_R_C', 'delta_slope_N2N3_C', 'SP_FFT_all_C', 'SP_AMP_all_C', 'SP_CDENS_all_C', 'SP_CHIRP_all_C', 'SP_COUPL_MAG_all_C', 'SP_COUPL_OVERLAP_all_C', 'SP_DENS_all_C', 'SP_ISA_S_all_C', 'SP_AMP_slow_C', 'SP_CDENS_slow_C', 'SP_CHIRP_slow_C', 'SP_COUPL_MAG_slow_C', 'SP_COUPL_OVERLAP_slow_C', 'SP_DENS_slow_C', 'SP_ISA_S_slow_C', 'SP_AMP_fast_C', 'SP_CDENS_fast_C', 'SP_CHIRP_fast_C', 'SP_COUPL_MAG_fast_C', 'SP_COUPL_OVERLAP_fast_C', 'SP_DENS_fast_C', 'SP_ISA_S_fast_C', 'SO_SLOPE_NEG1_C', 'SO_SLOPE_POS1_C', 'SO_RATE_C', 'SO_SLOPE_C', 'SO_DUR_C', 'SO_P2P_C']

    exposures = list(df2.columns[1:])#[1:2]
    mediators = eeg_names
    outcomes = ["Teng3MSScore_V2"]#, "TrialsBTime_V2"]

    # define outcome model and mediator model
    #    outcome: dementia
    #    exposure: cytokine
    #    mediator: sleep pattern

    covariates_basic = ['Educ', 'Race_Black','Race_Asian', 'Race_Hispanic', 'Race_Other', "APOE4Count"]
    #covariates_v1 = ['Age_V1', 'BMI_V1', 'Smoking', 'Alcohol', 'Med_AntiInfl_V1']
    covariates_vs1 = ["Age_VS1", "BMI_VS1", "MH_HTN_VS1", "MH_Stroke_VS1", "MH_DB2_VS1", "Med_Benzo_VS1", "Med_Antidep_VS1", "Med_Zolpidem_VS1", "Med_Opiod_VS1", "Med_AntiInfl_VS1", 'AHI4']

    # deal with missing value
    # we can use the simplest approach: exclude any subjects with any missing value
    cols = ['ID'] + exposures + mediators + outcomes + covariates_basic + covariates_vs1# + covariates_v1
    df = df[cols]
    print(df.shape)
    df = df.dropna(ignore_index=True)
    print(df.shape)

    df['Educ2'] = (df.Educ>=5).astype(int)
    df['Race_Other2'] = ((df.Race_Hispanic+df.Race_Other)>0).astype(int)
    df['BMI_VS1_cat'] = df.BMI_VS1.apply(lambda x:get_bmi_cat(x))
    df['MH'] = ((df.MH_HTN_VS1 + df.MH_DB2_VS1 + df.MH_Stroke_VS1)>0).astype(int)
    df['Med_Sleep_VS1'] = ((df.Med_Benzo_VS1 + df.Med_Antidep_VS1 + df.Med_Zolpidem_VS1)>0).astype(int)
    df['Med_AntiInfl_VS1_2'] = ((df.Med_Opiod_VS1 + df.Med_AntiInfl_VS1)>0).astype(int)
    df['AHI_VS1_cat'] = df.AHI_VS1.apply(lambda x:get_ahi_cat(x))

    caliper = {
    'Educ2':0
    'Race_Black':0,
    'Race_Asian':0,
    'Race_Other2':0,
    "APOE4Count":1,
    'Age_VS1':5,
    'BMI_VS1_cat':0,
    'MH':0,
    'Med_Sleep_VS1':0,
    'Med_AntiInfl_VS1_2':0,
    'AHI_VS1_cat':0,
    }
    


if __name__=='__main__':
    main()

