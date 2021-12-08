#!/usr/bin/env python
# coding: utf-8

#Venn diagrams
import pandas as pd
ibd_data = pd.read_pickle("datasetspandas/dutch_export_transposed_ibd_binary.pkl")
ibd_data.rename(columns={'Unnamed: 0': 'sample_id'}, inplace=True)
ibd_data.index.name = 'sample_id'

#Fetch cohort ID
covariates = pd.read_excel("Covariates/CDvsUC analysis covariates.xlsx")
covariates = covariates.set_index('sample_id')
del covariates['Plate_id']
del covariates['Age']
del covariates['Gender']

#Define groups
df_CD = covariates.loc[covariates['IBD type'] == 1]
df_UC = covariates.loc[covariates['IBD type'] == 0]
df_IBD = covariates.merge(ibd_data, on=["sample_id"])
cd_data = df_CD.merge(ibd_data, on=["sample_id"])
del cd_data['IBD type']
uc_data = df_UC.merge(ibd_data, on=["sample_id"])
del uc_data['IBD type']
del df_IBD['IBD type']

#CD feature selection (5% threshold)
def Prevalence(Series):
    Prev = Series.sum()/Series.count()
    return(Prev)
Prev_thresh = 0.05
print("Computing antibody prevalence")
Prevalences = df_IBD.apply(Prevalence, axis=0)
Vector_keep = Prevalences > Prev_thresh #True/False vector of the features with prevalence > 5%
df_IBD_f = df_IBD[Vector_keep.index[Vector_keep]] #Removal of columns that do not pass the prevalence threshold

#UC feature selection (5% threshold)
def Prevalence(Series):
    Prev = Series.sum()/Series.count()
    return(Prev)
Prev_thresh = 0.05
print("Computing antibody prevalence")
Prevalences = uc_data.apply(Prevalence, axis=0)
Vector_keep = Prevalences > Prev_thresh #True/False vector of the features with prevalence > 5%
uc_data_f = uc_data[Vector_keep.index[Vector_keep]] #Removal of columns that do not pass the prevalence threshold

cd_unique = cd_data_f.T
uc_unique = uc_data_f.T
cd_unique.to_excel("datasetspandas/cd_unique_features.xlsx")
uc_unique.to_excel("datasetspandas/uc_unique_features.xlsx")

df_IBD_f
ibd_unique = df_IBD_f.T
ibd_unique.to_excel("datasetspandas/ibd_unique_features.xlsx")

##Hereafter, random peptide selection exclusion and >95% exclusion was performed manually to compute overlapping numbers

#Figures
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
fig = plt.figure()

#Custom it
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set_style("white")
c = venn3_unweighted(subsets=(281, 243, "", "", 227, 122, 1974), set_labels = ('CD', 'UC', 'IBD'), alpha=1)
 
#Add title and annotation
c.get_patch_by_id('100')
c.get_patch_by_id('010')
c.get_patch_by_id('001')
plt.title("Prevalence of antibody probes (≥5%)", fontsize=18)
plt.savefig("datasetspandas/VennDiagramPrevalence-updatedsep2021.jpg", bbox_inches='tight', dpi=600)

from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#Custom it
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set_style("white")
c2 = venn2_unweighted(subsets = ('205', '104', '64'), set_labels = ('CD', 'UC'), alpha=0.5)
 
#Add title and annotation
c2.get_patch_by_id('10').set_color('deepskyblue')
c2.get_patch_by_id('01').set_color('lightblue')
c2.get_patch_by_id('11').set_color('darkblue')
c2.get_label_by_id('10').set_fontsize(18)
c2.get_label_by_id('01').set_fontsize(18)
c2.get_label_by_id('11').set_fontsize(18)
plt.title("Differentially abundant peptides", fontstyle="italic")
plt.suptitle("Total: 373 peptides (FDR≤0.05)", fontstyle="italic", x=0.52, y=0.78, fontsize=10)
plt.savefig("datasetspandas/VennDiagramCC-analysis-total_updatedaug2021.jpg", bbox_inches='tight', dpi=600)
