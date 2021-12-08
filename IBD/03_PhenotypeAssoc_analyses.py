#!/usr/bin/env python
# coding: utf-8


#Phenotype association analyses (IBD cohort-specific, e.g. surgical history, disease location, ASCA positivity, etc)
#Load Pandas and load antibody data
import pandas as pd
ibd_data = pd.read_pickle("datasetspandas/dutch_export_transposed_ibd_binary.pkl")
ibd_data

#Selection of prevalent antibodies
def Prevalence(Series):
    Prev = Series.sum()/Series.count()
    return(Prev)
Prev_thresh = 0.05
print("Computing antibody prevalence")
Prevalences = ibd_data.apply(Prevalence, axis=0)
Vector_keep = Prevalences > Prev_thresh #True/False vector of the features with prevalence > 5%
ibd_data_f = ibd_data[Vector_keep.index[Vector_keep]] #Removal of columns that do not pass the prevalence threshold

#Rename ID column to allow merging with covariates file
ibd_data_f = ibd_data_f.reset_index()
ibd_data_f.rename(columns={'Unnamed: 0': 'sample_id'}, inplace=True)
ibd_data_f

covariates = pd.read_excel("datasetspandas/Montreal-B1vsB2-CDcovariates.xlsx")

#Ensure numericity of data input and set indices similar for both dataframes to allow merging for input in logistic regression
ibd_data_f.apply(pd.to_numeric, errors='ignore')
covariates.apply(pd.to_numeric, errors='ignore')
covariates = covariates.set_index('sample_id')
ibd_data_f = ibd_data_f.set_index('sample_id')

#Load modules for running logistic regression
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

#Define function for logistic regression
def Do_logistic(ibd_data_f, Covariates, MontrealBehaviorCD):
   #Put covariates, phenotype and Y together and perform logistic regression 
    Name_probe = ibd_data_f.name
    #Putting Enrichment and regressors together
    ibd_data_f.name = "Enrichment"
    Input = pd.concat([ibd_data_f,covariates], axis=1)
    #Prepare formula
    Formula = "Enrichment ~ C(MontrealBehaviorCD) + Age + C(Gender)"
    Input["MontrealBehaviorCD"] = pd.to_numeric(Input["MontrealBehaviorCD"])
    Input["Gender"] = pd.to_numeric(Input["Gender"])
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

MontrealBehaviorCD = covariates['MontrealBehaviorCD']
Covariates = covariates[["Gender", "Age"]]
    
#df = Do_logistic(Data[e.g. "agilent_129495"], Covariates, IBD phenotype)
df = ibd_data_f.apply(Do_logistic, args=(Covariates, MontrealBehaviorCD),result_type="expand", axis=0)
df = df.T
df.to_csv("datasetspandas/resultsMontrealB1vsB2.csv")

#View results of logistic regression
results = pd.read_csv("datasetspandas/resultsMontrealB1vsB2.csv")

#Merge antibody and covariates df's to calculate antibody percentages per target group
datamerge = covariates.merge(ibd_data_f, on=["sample_id"])
#Remove unwanted covariates
del datamerge['Age']
del datamerge['Gender']
del datamerge['Plate_id']
datamerge

#Calculate percentages per target group and store them
cohort_count = datamerge.groupby(['MontrealBehaviorCD']).count()
cohort_sum = datamerge.groupby(['MontrealBehaviorCD']).sum()
cohort_percentages = (cohort_sum / cohort_count) * 100
cohort_percentages = cohort_percentages.T
cohort_percentages = cohort_percentages.reset_index()
cohort_percentages.to_excel("datasetspandas/cohort_percentages_B1vsB2.xlsx")

#Load percentages and rename columns to target group names
percentages = pd.read_excel("datasetspandas/cohort_percentages_B1vsB2.xlsx")
percentages = percentages.rename({0: '%MontrealB1'}, axis=1)
percentages = percentages.rename({1: '%MontrealB2'}, axis=1)

#Rename antibody peptide column to allow merging with results dataframe
del percentages['Unnamed: 0']
percentages = percentages.rename({'index': 'Unnamed: 0'}, axis=1)

#Merge percentages with results dataframe
resultswithpercentages = results.merge(percentages, on=["Unnamed: 0"])

#Load information on peptides of antibody dataframe
information = pd.read_csv("datasetspandas/df_info_AT.csv")

#Merge peptide info with results+percentages dataframe
results_proteininfo = resultswithpercentages.merge(information, on=["Unnamed: 0"])

#Remove unwanted variables
del results_proteininfo[{'aa_seq', 'len_seq', 'pos', 'is_pos_cntrl', 'is_neg_cntrl', 'is_phage', 'is_influenza', 'is_allergens', 'is_genome_editing', 'IEDB_DOIDs', 'IEDB_comments',
'IEDB_organism_name', 'IEDB_parent_species', 'is_rand_cntrl', 'is_VFDB', 'is_patho_strain', 'is_IgA_coated_strain', 'is_probio_strain', 'bac_src', 'is_gut_microbiome', 'Unnamed: 0'}]

results_proteininfo.to_excel("datasetspandas/MontrealBehaviorCDB1vsB2results.xlsx")

##Exclusion of random peptides due to identical sequences done after this step (see suppl tables)

#Plotting
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sns

colors_list = #specify colors here

sns.set(rc={'figure.figsize':(5.5,5.5)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
sns.set_style("white")
Plot = sns.scatterplot(x="%MontrealB1", y="%MontrealB2", hue="P-value ≤ 0.05", s=50, palette=colors_list, data=df)
plt.xlabel("% in non-stricturing, non-penetrating CD", size=14, labelpad=10)
plt.ylabel("% in stricturing CD", size=14, labelpad=10)
plt.title("Montreal B1 vs. Montreal B2 (CD)", size=14, y=1.03, fontweight='bold')

Plot.figure.savefig("Results/MontrealBehaviorCDB1vsB2.jpg", dpi=600, bbox_inches='tight')
