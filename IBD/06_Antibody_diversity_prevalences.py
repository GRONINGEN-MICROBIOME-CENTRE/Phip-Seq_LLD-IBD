#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

ibd_data = pd.read_pickle("datasetspandas/dutch_export_transposed_ibd_binary.pkl")

#Selection of prevalent features
def Prevalence(Series):
    Prev = Series.sum()/Series.count()
    return(Prev)
Prev_thresh = 0.05
print("Computing antibody prevalence")
Prevalences = ibd_data.apply(Prevalence, axis=0)
Vector_keep = Prevalences > Prev_thresh #True/False vector of the features with prevalence > 5% < 95%
ibd_data_f = ibd_data[Vector_keep.index[Vector_keep]] #Removal of columns that do not pass the prevalence threshold

ibd_data = ibd_data.reset_index()
ibd_data.rename(columns={'sample_id': 'Sample_id_WIS'}, inplace=True)

covariates = pd.read_excel("datasetspandas/DiseaseDurationcovariates.xlsx")
ibd_data.apply(pd.to_numeric, errors='ignore')
covariates.apply(pd.to_numeric, errors='ignore')

ibd_data = ibd_data.set_index('Sample_id_WIS')
ibd_data["sum_antibodies"] = ibd_data.sum(axis=1)

#Ensure to include: age, sex, plate ID, IBD type and disease duration
covariates = pd.read_excel("datasetspandas/Metadata IBD cohort_Arno_updated_aug2021.xlsx")
datamerge = covariates.merge(ibd_data, on=["Sample_id_WIS"])
dysbiosis = pd.read_csv("datasetspandas/Dysbiosis.score.txt", sep=" ")
dysbiosis.rename(columns={'ID': 'Participant ID'}, inplace=True)
datamergedysbiosis = dysbiosis.merge(datamerge, on=["Participant ID"])

#Dysbiosis vs antibody diversity
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(7,7)})
sns.set_style("white")
plot = sns.jointplot(data=datamergedysbiosis,x="sum_antibodies", y="Dysbiosis", color='darkblue', kind='reg', space=0)
plot.plot_marginals(sns.kdeplot, fill=True)
plot.set_axis_labels('No. of different antibodies', 'Dysbiosis score', fontsize=12)
plot.ax_marg_x.set_xlim(150, 2400)
plot.ax_marg_y.set_ylim(35, 65)

plot.savefig("datasetspandas/Microbiome-Dysbiosis-AntibodyDiversity-jointplot.jpg", dpi=600, bbox_inches='tight')

#Antibody prevalences
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

ibd_data["antibody_totals"] = ibd_data.sum(axis=0)
del ibd_data["sum_antibodies"]
del ibd_data["antibody_totals"]
ibd_data_totals = ibd_data.sum(axis=0)
ibd_data_totals_sorted = ibd_data_totals.sort_values(ascending=True)

df = pd.DataFrame({'prevalence':ibd_data_totals_sorted})
df = df.sort_values(by=['prevalence'])
numbers = ibd_data_totals_sorted
df["numbers"] = numbers

hue_order = list(range(1, 498))

#Rainbow plot of prevalences
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(7,7)})
sns.set_style("white")
sns.displot(x="prevalence", hue="prevalence", hue_order=hue_order, palette="Spectral", legend=False, data=df, bins=100)
plt.xlim(0, 501)
plt.xlabel('Number of patients', fontsize=12)
plt.ylabel('Number of present antibodies', fontsize=12)
plt.yscale('log')
plt.axhline(y=2368, xmin=0.05, linewidth=1, ls='--', color='black')
plt.axvline(x=26, ymin=0.65, ymax=0.671, linewidth=1, ls='-', color='black')
plt.axvline(x=499, ymin=0.65, ymax=0.671, linewidth=1, ls='-', color='black')
x1 = 275
y1 = 3000
plt.text(x1, y1, "2368 antibody-bound peptides in 5-95% of patients", ha='center', va='bottom', fontsize=9, fontstyle='italic')
plt.savefig("datasetspandas/PublicPrivateantibodiesIBDcohort-updated.jpg", dpi=600, bbox_inches='tight')

#Antibody diversity for CD and UC
sum_antibodies = ibd_data.sum(axis=1)
df2 = pd.DataFrame({'antibodies':sum_antibodies})
df2 = df2.sort_values(by=['antibodies'])
covariates = pd.read_excel("datasetspandas/DiseaseDurationcovariates.xlsx")
df2 = df2.reset_index()
df2.rename(columns={'Unnamed: 0': 'sample_id'}, inplace=True)
datamerge = covariates.merge(df2, on=["sample_id"])

df_CD = datamerge.loc[datamerge["IBDtype"] == 1]
df_UC = datamerge.loc[datamerge["IBDtype"] == 0]

from scipy.stats import mannwhitneyu
df_CD = datamerge.loc[datamerge["IBDtype"] == 1]
df_UC = datamerge.loc[datamerge["IBDtype"] == 0]
stat, p = mannwhitneyu(df_CD['antibodies'], df_UC['antibodies'])
print('Statistics=%.3f, p=%.3f' % (stat, p))
alpha = 0.05
if p > alpha:
    print('nonsignificant')
else:
    print('significant')

colors_list = ["#4B0082", "#DC143C"]

#Histograms for CD/UC
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(7,7)})
sns.set_style("white")
plot = sns.displot(data=datamerge, x="antibodies", hue='IBDtype', palette=colors_list, kind="hist", kde=True, legend=False, bins=70)
plot.set_axis_labels('No. of different antibodies', 'Frequency', fontsize=12)
plt.legend(title='Diagnosis', labels=['CD', 'UC'])
plt.xlim(0,3000)
plt.axvline(x=1015, linewidth=1, ls='--', color='black')
plt.axvline(x=720.5, linewidth=1, ls='--', color='grey')
plt.axvline(x=1362, linewidth=1, ls='--', color='grey')
x1, x2, x3 = 1140, 760, 1400
y1, y2, y3 = 13.7, 13.7, 13.7
plt.text(x1, y1, "Median", ha='center', va='bottom', rotation=35, fontsize=9)
plt.text(x2, y2, "Q1", ha='center', va='bottom', rotation=35, fontsize=9)
plt.text(x3, y3, "Q3", ha='center', va='bottom', rotation=35, fontsize=9)

plot.savefig("datasetspandas/Antibody_diversity_bydiagnosis.jpg", dpi=600, bbox_inches='tight')
