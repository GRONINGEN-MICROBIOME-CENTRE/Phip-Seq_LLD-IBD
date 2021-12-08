#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sns

#Prepare data, making it binary (0/1)
df = pd.read_pickle("datasetspandas/dutch_export_transposed.pkl")
ibd = df.loc[df['Unnamed: 0'].str.contains('33_')]
ibd.apply(lambda row: row.astype(str).str.contains('-1').any(), axis=1)
ibd2 = ibd.fillna(0)
ibd2 = ibd2.set_index('Unnamed: 0')
ibd3 = ibd2.apply(pd.to_numeric)
ibd3 = ibd3.where(ibd3 <= 0, 1)
ibd3 = ibd3.where(ibd3 >= 0, 0)

#Check if there any remaining NaNs (-1)
ibd3[ibd3.eq(-1).any(1)]

#Save
ibd3.to_pickle("datasetspandas/dutch_export_transposed_ibd_binary.pkl")

#PCA
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

#Read data
data = pd.read_pickle("datasetspandas/dutch_export_transposed_ibd_binary.pkl")

#Selection of prevalent features, removal of antibodies that are present in less than 5% of samples
def Prevalence(Series):
    Prev = Series.sum()/Series.count()
    return(Prev)
Prev_thresh = 0.05
print("Computing antibody prevalence")
Prevalences = data.apply(Prevalence, axis=0)
Vector_keep = Prevalences > Prev_thresh #True/False vector of the features with prevalence > 5%
data_f = data[Vector_keep.index[Vector_keep]] #Removal of columns that do not pass the prevalence threshold

#Perform PCA
print("Computing PCA")
pca = PCA(n_components=10) #specify no of principal components
pca.fit(data_f)
Egenvectors = pca.components_
Explained = pca.explained_variance_ratio_ * 100 #Percentage variance explained
Explained = [round(Explained[0],1)], round(Explained[1],1)
PCs = pca.transform(data_f) #Compute PCs using the eigenvectors calculated in the PCA function
Eigenvectors = Egenvectors.T
Eigenvectors = pd.DataFrame(Eigenvectors,columns=["Load_PC1", "Load_PC2", "Load_PC3", "Load_PC4", "Load_PC5", "Load_PC6", "Load_PC7", "Load_PC8", "Load_PC9", "Load_PC10"])
Eigenvectors["Probe"] = data_f.columns

#import sample info
metadata = pd.read_excel("datasetspandas/Metadata IBD cohort_Arno_updated_aug2021.xlsx")
metadata.set_index('Sample_id_WIS')

#Creating a pandas dataframe with the PCs, IDs and cohort info
PCs = pd.DataFrame(PCs,columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"])
PCs.index = data_f.index
data_f.index.names = ["Sample_id_WIS"]
PCs_merged = PCs.merge(metadata, on=["Sample_id_WIS"])

#In-between: examining loading factors of PCs
loadingfactors = pd.read_csv("datasetspandas/Loadings_PCIBD.csv")
loadingfactors.sort_values(by='Load_PC1', ascending=False)

proteindata = pd.read_csv("datasetspandas/df_info_AT.csv")
proteindata = proteindata.rename({'Unnamed: 0': 'Probe'}, axis=1)
loadingfactorswithproteins = loadingfactors.merge(proteindata, on=["Probe"])
del loadingfactorswithproteins[{'aa_seq', 'len_seq', 'pos', 'is_pos_cntrl', 'is_neg_cntrl', 'is_phage', 'is_influenza', 'is_allergens', 'is_genome_editing', 'IEDB_DOIDs', 'IEDB_comments',
'IEDB_organism_name', 'IEDB_parent_species', 'is_rand_cntrl', 'is_VFDB', 'is_patho_strain', 'is_IgA_coated_strain', 'is_probio_strain', 'bac_src', 'is_gut_microbiome', 'Unnamed: 0'}]
loadingfactorswithproteins.to_excel("datasetspandas/loadingfactorsIBDcohort.xlsx")

#Visualization of loading factors
sns.set(rc={'figure.figsize':(11.7,8.27)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
sns.set_style("white")
Plot = sns.scatterplot(x="Load_PC1", y="Load_PC2", data=loadingfactors)
#Add % of variance explained by the PC in each axis
Plot.set(xlabel="Load_PC1({N}%)".format(N=str(Explained[0])) , ylabel ="Load_PC2({N}%)".format(N=str(Explained[1])))
Plot.figure.savefig("datasetspandas/PCAIBD-loadingfactors.png")

#Clustering algorithm (KMeans). As visual inspection of PCA seems to indicate two major clusters, set k to 2
print("2-means clustering algorithm")
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=2, random_state=0).fit(data_f)
Clusters = kmeans.labels_
#Adding clusters to dataframe
Clusters = pd.DataFrame(Clusters,columns=["Cluster"])
PCs_merged = pd.concat([PCs_merged, Clusters],axis=1)

#Removing MC and DR from the variable 'IBD type' as they are meaningless
PCs_merged = PCs_merged[(PCs_merged['IBD type'] != 'MC') & (PCs_merged['IBD type'] != 'DR')]

#Removing IBDU from the variable 'IBD type' if required
PCs_merged = PCs_merged[(PCs_merged['IBD type'] != 'IBDU')]

#Coloring PCA plot based on probe prevalence / number of different antibodies per IBD patient
data_f["Probe_prevalence"] = data_f.sum(axis=1)
Sum_antibodies = data_f["Probe_prevalence"]
PCs_merged = PCs_merged.merge(Sum_antibodies, on=["Sample_id_WIS"], how='left')

#Plotting

#Plotting PCA colored by IBD subtype (without clustering algorithm)
colors_list = ["#DC143C", "#4B0082"]
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.scatterplot(x="PC1", y="PC2", hue="IBD type", palette=colors_list, data=PCs_merged)
#Add % of variance explained by the PC in each axis
Plot.set(xlabel="PC1({N}%)".format(N=str(Explained[0])) , ylabel ="PC2({N}%)".format(N=str(Explained[1])))
Plot.legend(title='Diagnosis')
Plot.figure.savefig("datasetspandas/PCAIBD-diagnosis_updatedaug2021.png", dpi=600, bbox_inches='tight')

#Plotting PCA colored by antibody diversity
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
sns.set_style("white")
Plot = sns.scatterplot(x="PC1", y="PC2", hue="Probe_prevalence", palette="rocket_r", data=PCs_merged)
#Add % of variance explained by the PC in each axis
Plot.set(xlabel="PC1({N}%)".format(N=str(Explained[0])) , ylabel ="PC2({N}%)".format(N=str(Explained[1])))
Plot.legend(title='Enriched antibodies', loc='upper right')
Plot.figure.savefig("datasetspandas/PCAIBD-Probe_prevalence_updatedaug2021.png", dpi=600, bbox_inches='tight')

#Plotting PCA colored by K-means cluster
colors_list = ['red', 'green']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
sns.set_style("white")
Plot = sns.scatterplot(x="PC1", y="PC2", hue="Cluster", palette=colors_list, data=PCs_merged)
#Add % of variance explained by the PC in each axis
Plot.set(xlabel="PC1({N}%)".format(N=str(Explained[0])) , ylabel = "PC2({N}%)".format(N=str(Explained[1])))
Plot.figure.savefig("datasetspandas/PCA_cluster_IBD_updatedaug2021.jpg", dpi=600, bbox_inches='tight')

#Test differences in PC1/2 between diagnosis and cluster identities
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

df_CD = PCs_merged.loc[PCs_merged["IBD type"] == 'CD']
df_UC = PCs_merged.loc[PCs_merged["IBD type"] == 'UC']
stat, p = mannwhitneyu(df_CD['PC1'], df_UC['PC1'])
print('Statistics=%.3f, p=%.3f' % (stat, p))
alpha = 0.05
if p > alpha:
    print('non-significant')
else:
    print('significant')

df_clusterzero = PCs_merged.loc[PCs_merged["Cluster"] == 0]
df_clusterone = PCs_merged.loc[PCs_merged["Cluster"] == 1]
stat, p = mannwhitneyu(df_clusterzero['PC2'], df_clusterone['PC2'])
print('Statistics=%.3f, p=%.3f' % (stat, p))
alpha = 0.05
if p > alpha:
    print('non-significant')
else:
    print('significant')

#Plotting PCA colored by plate ID (without clustering algorithm)
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.scatterplot(x="PC1", y="PC2", hue="Plate_id", data=PCs_merged)
#Add % of variance explained by the PC in each axis
Plot.set(xlabel="PC1({N}%)".format(N=str(Explained[0])) , ylabel ="PC2({N}%)".format(N=str(Explained[1])))
Plot.legend(title='Plate ID')
Plot.figure.savefig("datasetspandas/PCA-PlateID_updatedaug2021.png", dpi=600, bbox_inches='tight')

#Plotting PCA colored by serological CMV test results (without clustering algorithm)
PCs_merged["CMV-IgG-positivity"].replace({0: "Negative", 1: "Positive"}, inplace=True)
colors_list = ['darkturquoise', 'crimson']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.scatterplot(x="PC1", y="PC2", hue="CMV-IgG-positivity", palette=colors_list, data=PCs_merged)
#Add % of variance explained by the PC in each axis
Plot.set(xlabel="PC1({N}%)".format(N=str(Explained[0])) , ylabel ="PC2({N}%)".format(N=str(Explained[1])))
Plot.legend(title='Serological CMV status')
Plot.figure.savefig("datasetspandas/PCA-CMVstatus_updatedaug2021.png", dpi=600, bbox_inches='tight')

#Label PCA based on dysbiosis scores from 137 IBD pt
dysbiosis = pd.read_csv("datasetspandas/Dysbiosis.score.txt", sep=" ")
dysbiosis = dysbiosis.rename({'ID': 'Participant ID'}, axis=1)
dysbiosis = dysbiosis.set_index("Participant ID")
PCs_merged = PCs_merged.set_index("Participant ID")
PCs_merged_dysbiosis = PCs_merged.merge(dysbiosis, on=["Participant ID"])

#Plotting PCA colored by dysbiosis score in subset of 137 IBD patients (of whom MGS data were available within 1 year of sampling)
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.scatterplot(x="PC1", y="PC2", hue="Dysbiosis", palette="rocket_r", data=PCs_merged_dysbiosis)
#Add % of variance explained by the PC in each axis
Plot.set(xlabel="PC1", ylabel="PC2")
Plot.legend(title='Dysbiosis')
Plot.figure.savefig("datasetspandas/PCAIBD-dysbiosisscore_sep2021.png", dpi=600, bbox_inches='tight')

#Examine associations between PCs and factors using linear regression
from sklearn.linear_model import LinearRegression

InLR2 = PCs_merged[['PC2', 'Probe_prevalence']]
Y = InLR2.iloc[:, 0].values.reshape(-1, 1)
X = InLR2.iloc[:, 1].values.reshape(-1, 1)
linear_regressor = LinearRegression()
linear_regressor.fit(X, Y)
Y_pred = linear_regressor.predict(X)
print('intercept:', linear_regressor.intercept_)
print('slope:', linear_regressor.coef_)

plt.scatter(X, Y)
plt.plot(X, Y_pred, color='red')
plt.ylabel('PC1')
plt.xlabel('Plate_id')
plt.legend(['PC1'])
plt.title('Plate-ID vs. PC1')
plt.grid()
plt.savefig("datasetspandas/LRplateidPC1.png")
plt.show()

LinReg = PCs_merged_dysbiosis[['PC10', 'Dysbiosis']]
Y = LinReg.iloc[:, 0].values.reshape(-1, 1)
X = LinReg.iloc[:, 1].values.reshape(-1, 1)

import statsmodels.api as sm
from scipy import stats
import statsmodels.formula.api as smf
linear_model=sm.OLS(Y,X)
result=linear_model.fit()
print(result.summary())

# standardizing dataframe
df_z = LinReg.select_dtypes(include=[np.number]).dropna().apply(stats.zscore)

# fitting regression
formula = 'PC2 ~ Plate_id'
result = smf.ols(formula, data=df_z).fit()

# checking results
result.summary()

#In-between: checking relationships between PCs and the influence of plate ID's
#Plotting plate ID's against the first two PCs
sns.set(rc={'figure.figsize':(11.7,8.27)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
sns.set_style("white")
PCPlate = sns.scatterplot(x="Plate_id", y="PC2", data=PCs_merged)
PCPlate.set(xlabel="Plate ID", ylabel ="PC2")
PCPlate.figure.savefig("datasetspandas/PCAIBD_plateid_vs_PC2.png")

#Tree map visualisation of antigen library composition
import squarify

#Put numbers in df
df = pd.DataFrame({'oligo_number':[122551,24510,24500,22050,14700,24164,11525,40000,30000,30000], 'group':["Bacterial genes", "Bacterial strains", "Pathogenic bacteria", "Antibody-coated bacteria", "Probiotic bacteria", "Virulence factors", "Controls", "Phages", "Allergens", "Immune epitopes"]})

#Define color palette
pal = sns.color_palette('tab10', 10)
pal = pal.as_hex()
print(pal)

sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
sns.set(rc={'figure.figsize':(10,7)})
colors = ['#1f77b4', '#17becf', '#d62728', '#9467bd', '#2ca02c', '#ff7f0e', '#7f7f7f', '#e377c2','#bcbd22', '#8c564b']
plot = squarify.plot(sizes=df['oligo_number'], label=df['group'], color=colors, text_kwargs={'size':12, 'weight': 'bold'}, alpha=.8)
x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = 21, 21, 58, 87, 53.5, 82.5, 56, 56, 85, 85
y1, y2, y3, y4, y5, y6, y7, y8, y9, y10 = 37, 87, 7, 7, 29, 29, 45.4, 73.5, 81.5, 53
plt.text(x1, y1, "122,551 oligos (~36%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x2, y2, "24,510 oligos (~7%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x3, y3, "24,500 oligos (~7%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x4, y4, "22,050 oligos (~6%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x5, y5, "14,700 oligos (~4%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x6, y6, "24,164 oligos (~7%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x7, y7, "11,525 oligos (~3%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x8, y8, "40,000 oligos (~12%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x9, y9, "30,000 oligos (~9%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.text(x10, y10, "30,000 oligos (~9%)", ha='center', va='bottom', fontsize=8, fontstyle='italic')
plt.title("Antigen library content", size=14, fontweight='bold', y=1.02)
plt.axis('off')
plot.figure.savefig("datasetspandas/TreeMapLibraryplot_nov2021.jpg", dpi=600, bbox_inches='tight')
