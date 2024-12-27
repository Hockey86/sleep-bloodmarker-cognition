import numpy as np
import pandas as pd
from pingouin import partial_corr
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12.5})
import matplotlib.cm as cm


df = pd.read_csv('dataset_preprocessed.csv')
       
cols_cov = ['Educ', 'Race_Black', 'Race_Asian', 'Race_Hispanic',
'Race_Other', 'APOE4Count', 'Age_VS1', 'BMI_VS1', 'MH_HTN_VS1',
'MH_Stroke_VS1', 'MH_DB2_VS1', 'Med_Benzo_VS1', 'Med_Antidep_VS1',
'Med_Zolpidem_VS1', 'Med_Opiod_VS1', 'Med_AntiInfl_VS1', 'AHI4']
cols_biomarker = ['Adiponectin_VS1', 'Leptin_VS1', 'Glucose_VS1', 'Insulin_VS1', 'CRP_VS1', 'IFNg_VS1', 'IL6_VS1', 'TNFa_VS1', 'TNFaSR2_VS1']
cols_sleep = ['delta_dbs_N3_C', 'delta_rel_N3_C', 'alpha_dbs_N3_C',
'alpha_rel_N3_C', 'delta_dbs_N2_C', 'delta_rel_N2_C', 'theta_dbs_N1_C',
'theta_rel_N1_C', 'delta_dbs_R_C', 'delta_rel_R_C', 'theta_dbs_R_C',
'theta_rel_R_C', 'delta_slope_N2N3_C', 'SP_FFT_all_C', 'SP_AMP_all_C',
'SP_CDENS_all_C', 'SP_CHIRP_all_C', 'SP_COUPL_MAG_all_C',
'SP_COUPL_OVERLAP_all_C', 'SP_DENS_all_C', 'SP_ISA_S_all_C',
'SP_AMP_slow_C', 'SP_CDENS_slow_C', 'SP_CHIRP_slow_C',
'SP_COUPL_MAG_slow_C', 'SP_COUPL_OVERLAP_slow_C', 'SP_DENS_slow_C',
'SP_ISA_S_slow_C', 'SP_AMP_fast_C', 'SP_CDENS_fast_C',
'SP_CHIRP_fast_C', 'SP_COUPL_MAG_fast_C', 'SP_COUPL_OVERLAP_fast_C',
'SP_DENS_fast_C', 'SP_ISA_S_fast_C', 'SO_SLOPE_NEG1_C',
'SO_SLOPE_POS1_C', 'SO_RATE_C', 'SO_SLOPE_C', 'SO_DUR_C', 'SO_P2P_C',]
#For band powers, the effective number of comparisons was 6. For spindle and SOs, the effective number of comparisons was 2. 
col_cog = 'Teng3MSScore_V2'

corrs_cog_sleep = np.zeros(len(cols_sleep))
pvals_cog_sleep = np.zeros(len(cols_sleep))
for i, col in enumerate(cols_sleep):
    res = partial_corr(data=df, x=col_cog, y=col, covar=cols_cov, alternative='two-sided', method='spearman')
    corrs_cog_sleep[i] = res.r.iloc[0]
    pvals_cog_sleep[i] = res['p-val'].iloc[0]
sigs_cog_sleep = pvals_cog_sleep<0.05/np.array([6]*13+[2]*28)

corrs_cog_biomarker = np.zeros(len(cols_biomarker))
pvals_cog_biomarker = np.zeros(len(cols_biomarker))
for i, col in enumerate(cols_biomarker):
    res = partial_corr(data=df, x=col_cog, y=col, covar=cols_cov, alternative='two-sided', method='spearman')
    corrs_cog_biomarker[i] = res.r.iloc[0]
    pvals_cog_biomarker[i] = res['p-val'].iloc[0]
sigs_cog_biomarker = pvals_cog_biomarker<0.05/8

corrs_biomarker_sleep = np.zeros((len(cols_biomarker), len(cols_sleep)))
pvals_biomarker_sleep = np.zeros((len(cols_biomarker), len(cols_sleep)))
for i, col1 in enumerate(cols_biomarker):
    for j, col2 in enumerate(cols_sleep):
        res = partial_corr(data=df, x=col1, y=col2, covar=cols_cov, alternative='two-sided', method='spearman')
        corrs_biomarker_sleep[i,j] = res.r.iloc[0]
        pvals_biomarker_sleep[i,j] = res['p-val'].iloc[0]
sigs_biomarker_sleep = pvals_biomarker_sleep<0.05/8/np.array([6]*13+[2]*28)


ids = np.where(np.vstack([sigs_cog_sleep,sigs_biomarker_sleep]).any(axis=0))[0]
cols_sleep = [cols_sleep[i] for i in ids]
sigs_cog_sleep = sigs_cog_sleep[ids]
corrs_cog_sleep = corrs_cog_sleep[ids]
pvals_cog_sleep = pvals_cog_sleep[ids]
sigs_sleep_biomarker = sigs_biomarker_sleep[:,ids].T
corrs_sleep_biomarker = corrs_biomarker_sleep[:,ids].T
pvals_sleep_biomarker = pvals_biomarker_sleep[:,ids].T
cols_biomarker = ['Adiponectin', 'Leptin', 'Glucose', 'Insulin', 'CRP', r'IFN-$\gamma$', 'IL-6', r'TNF-$\alpha$', r'TNF-$\alpha$sRII']
assert cols_sleep==['delta_dbs_N3_C', 'theta_rel_N1_C', 'SP_FFT_all_C', 'SP_CDENS_all_C', 'SP_COUPL_OVERLAP_all_C', 'SP_DENS_all_C', 'SP_CDENS_slow_C', 'SP_COUPL_OVERLAP_slow_C', 'SP_DENS_slow_C', 'SP_CDENS_fast_C', 'SP_COUPL_OVERLAP_fast_C', 'SP_DENS_fast_C', 'SO_DUR_C']
cols_sleep = ['delta N3', 'relative theta N1', 'spindle frequency', 'SO-coupled spindle density', 'spindle-SO overlap', 'spindle density', 'SO-coupled slow spindle density', 'slow spindle-SO overlap', 'slow spindle density', 'SO-coupled fast spindle density', 'fast spindle-SO overlap', 'fast spindle density', 'SO duration']
print(corrs_sleep_biomarker.shape)
col_cog = '3MS\nat\nVisit 2'

print(f'sigs_sleep_biomarker = {sigs_sleep_biomarker}')
print(f'sigs_cog_biomarker = {sigs_cog_biomarker}')
print(f'sigs_cog_sleep = {sigs_cog_sleep}')

Ns, Nb = corrs_sleep_biomarker.shape
frame_color = (0.5,0.5,0.5)

def pval2radius(p):
    max_ = 0.45
    min_ = 0.1
    r = (min_-max_)/np.sqrt(0.05)*np.sqrt(p)+max_  # 0-->max, 0.05-->min
    return r


def effectsize2color(f, vmin, vmax, cmap_name):
    cmap = mpl.colormaps[cmap_name]
    c = cmap(np.clip((f-vmin)/(vmax-vmin), 0,1))
    return c

vmin = -0.1
vmax = 0.1
cmap_name = 'coolwarm'

plt.close()
fig = plt.figure(figsize=(7,6.5))
gs = fig.add_gridspec(2, 3, height_ratios=(1,Ns), width_ratios=(4,Nb,1))

ax_sleep_biomarker = fig.add_subplot(gs[1,1])
for i in range(Ns+1):
    ax_sleep_biomarker.plot([0,Nb], [Ns-i]*2, c=frame_color, lw=1)
    #if i<Ns:
    #    ax_sleep_biomarker.text(-0.2, Ns-i-0.5, cols_sleep[i], ha='right', va='center')
for i in range(Nb+1):
    ax_sleep_biomarker.plot([i,i], [0,Ns], c=frame_color, lw=1)
for i in range(Ns):
    for j in range(Nb):
        if pvals_sleep_biomarker[i,j]<0.05:
            r = pval2radius(pvals_sleep_biomarker[i,j])
            c = effectsize2color(corrs_sleep_biomarker[i,j], vmin, vmax, cmap_name)
            ax_sleep_biomarker.add_patch(plt.Circle((j+0.5,Ns-i-0.5), r, fc=c, ec='k'))
        if sigs_sleep_biomarker[i,j]:
            ax_sleep_biomarker.scatter([j+0.5], [Ns-i-0.5], s=15, color='k')
ax_sleep_biomarker.axis('off')

ax = fig.add_subplot(gs[1,0], sharey=ax_sleep_biomarker)
for i in range(Ns):
    ax.text(1, Ns-i-0.5, cols_sleep[i], ha='right', va='center')
ax.set_xlim(0,1)
ax.axis('off')

ax_cog_biomarker = fig.add_subplot(gs[0,1], sharex=ax_sleep_biomarker)
for i in range(Nb+1):
    ax_cog_biomarker.plot([i,i], [0,1], c=frame_color, lw=1)
    if i<Nb:
        ax_cog_biomarker.text(i+0.5, 1+0.2, cols_biomarker[i], ha='center', va='bottom', rotation=90)
for i in [0,1]:
    ax_cog_biomarker.plot([0,Nb], [i,i], c=frame_color, lw=1)
for i in range(Nb):
    if pvals_cog_biomarker[i]<0.05:
        r = pval2radius(pvals_cog_biomarker[i])
        c = effectsize2color(corrs_cog_biomarker[i], vmin, vmax, cmap_name)
        ax_cog_biomarker.add_patch(plt.Circle((i+0.5,0.5), r, fc=c, ec='k'))
    if sigs_cog_biomarker[i]:
        ax_cog_biomarker.scatter([i+0.5], [0.5], s=15, color='k')
ax_cog_biomarker.axis('off')

ax_cog_sleep = fig.add_subplot(gs[1,2], sharey=ax_sleep_biomarker)
for i in range(Ns+1):
    ax_cog_sleep.plot([0,1], [i,i], c=frame_color, lw=1)
for i in [0,1]:
    ax_cog_sleep.plot([i,i], [0,Ns], c=frame_color, lw=1)
for i in range(Ns):
    if pvals_cog_sleep[i]<0.05:
        r = pval2radius(pvals_cog_sleep[i])
        c = effectsize2color(corrs_cog_sleep[i], vmin, vmax, cmap_name)
        ax_cog_sleep.add_patch(plt.Circle((0.5,Ns-i-0.5), r, fc=c, ec='k'))
    if sigs_cog_sleep[i]:
        ax_cog_sleep.scatter([0.5], [Ns-i-0.5], s=15, color='k')
ax_cog_sleep.axis('off')

ax = fig.add_subplot(gs[0,0])
fig.colorbar(cm.ScalarMappable(
    norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True),
    cmap=cmap_name), ax=ax,
    orientation='vertical', shrink=3.6, aspect=10, anchor=(0.65,0), label='Correlation', location='left', fraction=1)
ax.axis('off')

ax = fig.add_subplot(gs[0,2])
ax.text(0,0,col_cog,ha='center',va='bottom')
ax.set_xlim(-1,1)
ax.set_ylim(0,1)
ax.axis('off')

plt.tight_layout()
plt.subplots_adjust(hspace=0, wspace=0)
#plt.show()
plt.savefig('correlations.pdf', bbox_inches='tight', pad_inches=0.01)#, dpi=300)

