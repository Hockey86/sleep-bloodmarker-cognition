import os, pickle
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from scipy.special import logit, expit
from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectFpr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ARDRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from skopt import BayesSearchCV
from skopt.space import Real, Categorical, Integer
from tqdm import tqdm
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import seaborn as sns
sns.set_style('ticks')


def filter_func(X, y):
    stats = []
    pvals = []
    for i in range(X.shape[1]):
        stat, pval = spearmanr(X[:,i],y)
        stats.append(stat)
        pvals.append(pval)
    return np.array(stats), np.array(pvals)
    
    
def get_perf(y, yp, nbt=0):
    rmse = []
    mae = []
    pearson_corr = []
    spearman_corr = []
    for i in tqdm(range(nbt+1)):
        btids = np.random.choice(len(y), len(y), replace=True)
        y2 = y[btids]
        yp2 = yp[btids]
        rmse.append(np.sqrt(np.mean((y2-yp2)**2)))
        mae.append(np.mean(np.abs(y2-yp2)))
        pearson_corr.append(pearsonr(y2,yp2)[0])
        spearman_corr.append(spearmanr(y2,yp2)[0])
        
    res = pd.Series()
    res['RMSE'] = rmse[0]
    res['MAE'] = mae[0]
    res['PearsonR'] = pearson_corr[0]
    res['SpearmanR'] = spearman_corr[0]
    res['N'] = len(y)
    if nbt>=10:
        res['RMSE_ub'] = np.percentile(rmse, 2.5)
        res['RMSE_lb'] = np.percentile(rmse, 97.5)
        res['MAE_ub'] = np.percentile(mae, 2.5)
        res['MAE_lb'] = np.percentile(mae, 97.5)
        res['PearsonR_ub'] = np.percentile(pearson_corr, 2.5)
        res['PearsonR_lb'] = np.percentile(pearson_corr, 97.5)
        res['SpearmanR_ub'] = np.percentile(spearman_corr, 2.5)
        res['SpearmanR_lb'] = np.percentile(spearman_corr, 97.5)
    
    return res
    
    
def fit_model(X, y, cv_K, model_type, random_state=None, n_jobs=1):
    y2 = logit(y/101)##
    inv_func = lambda x: np.clip(expit(x)*101, 0, 100)
    
    alpha = 1#0.4
    kf = KFold(n_splits=cv_K, random_state=random_state, shuffle=True)
    ypte_cv = np.zeros_like(y)
    cv_ids = np.zeros(len(y), dtype=int)
    models_cv = []
    for cvi, (trids, teids) in enumerate(tqdm(kf.split(y.reshape(-1,1)), total=cv_K)):
        Xtr = X[trids]; ytr = y2[trids]
        Xte = X[teids]
        
        scaler = StandardScaler()
        feat_selection = SelectFpr(filter_func, alpha=alpha)
        if model_type == 'lr':
            pred_model = ARDRegression(max_iter=300, tol=0.001, fit_intercept=True, verbose=False)
            hp = {}
            hp_dtype = {}
        elif model_type == 'rf':
            pred_model = RandomForestRegressor(random_state=random_state+cvi)
            hp = {
                'n_estimators':Integer(2,50),
                'max_depth': Integer(1,5),
                'min_samples_split': Integer(2,5),
                #'min_samples_leaf': Integer(2,5),
                'ccp_alpha': Real(1e-3,1e3, prior='log-uniform'),
                }
            hp_dtype = {
                'n_estimators':int,
                'max_depth': int,
                'min_samples_split': int,
                #'min_samples_leaf': int,
                'ccp_alpha': float,
                }
        else:
            raise NotImplementedError(f'model_type = {model_type}')

        model = Pipeline([
            ('scaler', scaler),
            ('feat_selection', feat_selection),
            ('pred_model', pred_model)])

        if len(hp)>0:
            model = BayesSearchCV(model,
                {'pred_model__'+k:v for k,v in hp.items()}, n_iter=50,
                scoring='neg_root_mean_squared_error',
                n_jobs=n_jobs, n_points=10, iid=True, refit=True,
                cv=cv_K, verbose=0, random_state=random_state,
                #error_score=-np.inf
                )
        
        model.fit(Xtr, ytr)
        ypte_cv[teids] = model.predict(Xte)
        models_cv.append(model)
        cv_ids[teids] = cvi
    ypte_cv = inv_func(ypte_cv)
    perf_cv = get_perf(y, ypte_cv, nbt=1000)
    
    # refit
    scaler = StandardScaler()
    feat_selection = SelectFpr(filter_func, alpha=alpha)
    if model_type == 'lr':
        pred_model = ARDRegression(max_iter=300, tol=0.001, fit_intercept=True, verbose=False)
    elif model_type == 'rf':
        hps = {}
        for hp_name in hp.keys():
            dt = hp_dtype[hp_name]
            hps[hp_name] = dt(np.median([x.best_params_['pred_model__'+hp_name] for x in models_cv]))
        pred_model = RandomForestRegressor(random_state=random_state+cvi, **hps)
    else:
        raise NotImplementedError(f'model_type = {model_type}')

    model = Pipeline([
        ('scaler', scaler),
        ('feat_selection', feat_selection),
        ('pred_model', pred_model)])

    model.fit(X, y2)
    yp = model.predict(X)
    yp = inv_func(yp)
    perf = get_perf(y, yp, nbt=1000)
    feat_support = model.named_steps['feat_selection'].get_support()
        
    if len(hp)>0:
        models_cv = [x.best_estimator_ for x in models_cv]
    return perf_cv, models_cv, ypte_cv, cv_ids, perf, model, yp, feat_support
    

def main():
    df = pd.read_csv('dataset_preprocessed.csv')
    
    output_dir = 'prediction_results'
    os.makedirs(output_dir, exist_ok=True)
    
    yname = "Teng3MSScore_V2"
    cov_names = list(df.columns)[52:]
    sleep_names = list(df.columns)[1:42]
    biomarker_names = list(df.columns)[42:51]
    
    random_state = 2024
    y = df[yname].values.astype(float)
    K = 10
    col_names = ['cov', 'biomarker', 'sleep', 'sleep+biomarker', 'cov+biomarker', 'cov+sleep', 'cov+sleep+biomarker']
    model_types = ['lr']#, 'rf']
    n_jobs = 14#
    
    for ni, name in enumerate(col_names):
        print(name)
        if name=='sleep+biomarker':
            cols = sleep_names+biomarker_names
        elif name=='sleep':
            cols = sleep_names
        elif name=='biomarker':
            cols = biomarker_names
        elif name=='cov+sleep+biomarker':
            cols = cov_names+sleep_names+biomarker_names
        elif name=='cov+sleep':
            cols = cov_names+sleep_names
        elif name=='cov+biomarker':
            cols = cov_names+biomarker_names
        elif name=='cov':
            cols = cov_names
            
        for model_type in model_types:
            print(model_type)
            output_path = os.path.join(output_dir, f'results_{name}_{model_type}.pickle')
            X = df[cols].values
            perf_cv, models_cv, ypte_cv, cv_ids, perf, model, yp, feat_support = fit_model(X, y, K, model_type, random_state=random_state+ni, n_jobs=n_jobs)
            print(perf_cv)
            cols_ = np.array(cols)[feat_support]
            df_ = pd.DataFrame(data={'Name':cols_, 'Value':model.named_steps['pred_model'].coef_.flatten()})
            df_ = df_[np.abs(df_.Value)>0].reset_index(drop=True)
            print(df_)
            with open(output_path, 'wb') as ff:
                pickle.dump({'sids':df.ID.values, 'y':y,
                'perf_cv':perf_cv, 'models_cv':models_cv,
                'ypte_cv':ypte_cv, 'cv_ids':cv_ids,
                'perf':perf, 'model':model, 'yp':yp,
                'names':cols, 'names_selected':cols_}, ff)
    
    model_type = 'lr'
    df_perf = []
    for ni, name in enumerate(col_names):
        path = os.path.join(output_dir, f'results_{name}_{model_type}.pickle')
        with open(path, 'rb') as ff:
            res = pickle.load(ff)
        ypte_cv = res['ypte_cv']
        #y = res['y']
        perf_cv = res['perf_cv']
        df_perf.append(perf_cv)
        
        plt.close()
        fig = plt.figure(figsize=(6,30/8))
        ax = fig.add_subplot(111)
        ax.scatter(y, ypte_cv, s=3, c='k')
        ax.plot([60,100], [60,100], c='r', ls='--')
        ax.set_xlabel('3MS at Visit 2')
        ax.set_ylabel('Predicted 3MS at Visit 2')
        ax.set_xlim(60,101)
        ax.set_ylim(75,101)
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'scatterplot_{name}_{model_type}.png'), bbox_inches='tight')
    
    df_perf = pd.DataFrame(df_perf)
    df_perf.insert(0, 'Name', col_names)
    df_perf.to_excel(os.path.join(output_dir, f'perf_{model_type}.xlsx'), index=False)


if __name__=='__main__':
    main()
    
