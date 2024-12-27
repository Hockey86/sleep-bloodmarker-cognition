import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


df = pd.read_csv('dataset_preprocessed.csv')
cols = list(df.columns)

sleep_cols1 = [x for x in cols if '_dbs_' in x] + ['delta_slope_N2N3_C']
sleep_cols2 = [x for x in cols if '_rel_' in x]
sleep_cols3 = ['SP_AMP_all_C', 'SP_CDENS_all_C', 'SP_DENS_all_C', 'SP_ISA_S_all_C', 'SP_COUPL_MAG_all_C', 'SP_COUPL_OVERLAP_all_C']
sleep_cols4 = [x for x in cols if 'SO_' in x] # --> should be 2 effective dimensions, duration and amplitude
biomarker_cols = cols[cols.index('Adiponectin_VS1'):cols.index('Insulin_VS1')+1]


for cols in [sleep_cols1, sleep_cols2, sleep_cols3, sleep_cols4, biomarker_cols]:
    X = df[cols].values
    X = X[~np.isnan(X).any(axis=1)]
    #X = (X-X.mean(axis=0))/X.std(axis=0)
    X2 = PCA(n_components=0.99).fit_transform(X)
    print(X.shape[1], X2.shape[1])

