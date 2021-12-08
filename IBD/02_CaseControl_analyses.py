#!/usr/bin/env python
# coding: utf-8

#Case-control analyses, comparing peptide frequencies between CD, UC and LL while adjusting for effects of age & gender

#Load Pandas and load antibody data
import pandas as pd
ibd_ll_data = pd.read_pickle("datasetspandas/dutch_export_transposed_ibdll_binary.pkl")

#Selection of prevalent features, removal of antibodies that are present in less than 5% of samples (and >95% of samples)
def Prevalence(Series):
    Prev = Series.sum()/Series.count()
    return(Prev)
Prev_thresh = 0.05 #0.95 to exclude peptides >95%
print("Computing antibody prevalence")
Prevalences = ibd_ll_data.apply(Prevalence, axis=0)
Vector_keep = Prevalences > Prev_thresh #True/False vector of the features with prevalence > 5%, convert to exclude > 95% peptides
ibd_ll_f = ibd_ll_data[Vector_keep.index[Vector_keep]] #Removal of columns that do not pass the prevalence threshold

#Include 323 extra probes that are >5% in LL but <5% in IBD
df_extraCC = pd.read_excel("datasetspandas/List_antibodies_to_add_CC.xlsx")

#Concatenate peptides to existing df
CCdf = pd.concat([ibd_ll_f, df_extraCC], axis=1)

#Rename ID column to allow merging with covariates file
CCdf = CCdf.reset_index()
CCdf.rename(columns={'Unnamed: 0': 'sample_id'}, inplace=True)

#Load covariates file
covariates = pd.read_excel("datasetspandas/CC_CD_covariates_March2021_matched.xlsx")

#Ensure numericity of data input and set indices similar for both dataframes to allow merging for input in logistic regression
CCdf.apply(pd.to_numeric, errors='ignore')
covariates.apply(pd.to_numeric, errors='ignore')
covariates = covariates.set_index('sample_id')
CCdf = CCdf.set_index('sample_id')

#Load necessary packages for running logistic regression
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

#Define function for logistic regression
def Do_logistic(CCdf, Covariates, Cohort):
   #Put covariates, phenotype and Y together and perform logistic regression 
    Name_probe = CCdf.name
    #Putting Enrichment and regressors together
    CCdf.name = "Enrichment"
    Input = pd.concat([CCdf,covariates], axis=1)
    #Prepare formula
    Formula = "Enrichment ~ C(Cohort) + Age + C(Sex)"
    Input["Cohort"] = pd.to_numeric(Input["Cohort"])
    Input["Sex"] = pd.to_numeric(Input["Sex"])
    Input["Age"] = pd.to_numeric(Input["Age"])
    #Remove samples with NA
    Input = Input.dropna(axis=0)
    #Call function
    log_reg = smf.logit(formula=Formula,data=Input).fit(maxiter=100,method="bfgs")
    #except: return(Name, "NA", "NA")
    Pvalue = log_reg.pvalues[1]
    Estimate = log_reg.params[1]
    Out = [Name_probe, Pvalue, Estimate]
    return(Out)

Cohort = covariates['Cohort']
Covariates = covariates[["Sex", "Age"]]
    
#df = Do_logistic(Data[e.g. "agilent_129495"], Covariates, IBD type)
df = CCdf.apply(Do_logistic, args=(Covariates, Cohort),result_type="expand", axis=0)
df = df.T
df.to_csv("datasetspandas/CDvsHC_fullresults_November2021_CCmatching.csv")

#View results of logistic regression
results = pd.read_csv("datasetspandas/CDvsHC_fullresults_November2021_CCmatching.csv")

#Merge antibody and covariates df's to calculate antibody percentages per target group
datamerge = covariates.merge(CCdf, on=["sample_id"])

#Remove unwanted covariates
del datamerge['Age']
del datamerge['Sex']

#Calculate percentages per target group and store them
cohort_count = datamerge.groupby(['Cohort']).count()
cohort_sum = datamerge.groupby(['Cohort']).sum()
cohort_percentages = (cohort_sum / cohort_count) * 100
cohort_percentages = cohort_percentages.T
cohort_percentages = cohort_percentages.reset_index()
cohort_percentages.to_excel("datasetspandas/cohort_percentages_CDvsHC_November2021_CCmatching.xlsx")

#Load percentages and rename columns to target group names
percentages = pd.read_excel("datasetspandas/cohort_percentages_CDvsHC_November2021_CCmatching.xlsx")
percentages = percentages.rename({0: '%LL'}, axis=1)
percentages = percentages.rename({1: '%CD'}, axis=1)

#Rename antibody probe column to allow merging with results dataframe
del percentages['Unnamed: 0']
percentages = percentages.rename({'index': 'Unnamed: 0'}, axis=1)

#Merge percentages with results dataframe
resultswithpercentages = results.merge(percentages, on=["Unnamed: 0"])

#Load information on proteins of antibody dataframe
information = pd.read_csv("datasetspandas/df_info_AT.csv")

#Merge protein info with results+percentages dataframe
results_proteininfo = resultswithpercentages.merge(information, on=["Unnamed: 0"])

#Remove unwanted variables
del results_proteininfo[{'aa_seq', 'len_seq', 'pos', 'is_pos_cntrl', 'is_neg_cntrl', 'is_phage', 'is_influenza', 'is_allergens', 'is_genome_editing', 'IEDB_DOIDs', 'IEDB_comments',
'IEDB_organism_name', 'IEDB_parent_species', 'is_rand_cntrl', 'is_VFDB', 'is_patho_strain', 'is_IgA_coated_strain', 'is_probio_strain', 'bac_src', 'is_gut_microbiome', 'Unnamed: 0'}]

#Store the final dataframe with results, percentages of target groups and protein annotation
results_proteininfo.to_excel("datasetspandas/CDvsHC_fullresults_November2021_CCmatching.xlsx")

#Plotting example
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sns
import pandas as pd

figuredata = pd.read_excel("datasetspandas/CDvsHC_fullresults_November2021.xlsx")

colors_list = ['salmon', 'dodgerblue']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.scatterplot(x="%LL", y="%CD", hue="Significant (FDR<0.05)", s=50, palette=colors_list, data=figuredata)
plt.xlabel("% in HC", size=16, labelpad=10)
plt.ylabel("% in CD", size=16, labelpad=10)
plt.title("CD vs. HC", size=16, y=1.03, fontweight='bold')
Plot.figure.savefig("datasetspandas/CDvsHC_March2021-updatednov2021.jpg")

