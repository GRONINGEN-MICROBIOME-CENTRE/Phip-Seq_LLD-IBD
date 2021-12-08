#!/usr/bin/env python
# coding: utf-8

#Prediction/classification part

#Load pandas and load antibody data of IBD and LL combined
import pandas as pd
ibd_ll_data = pd.read_pickle("datasetspandas/dutch_export_transposed_ibdll_binary.pkl")
ibd_ll_data

#Selection of prevalent features, removal of antibodies that are present in less than 5% of samples
def Prevalence(Series):
    Prev = Series.sum()/Series.count()
    return(Prev)
Prev_thresh = 0.05
print("Computing antibody prevalence")
Prevalences = ibd_ll_data.apply(Prevalence, axis=0)
Vector_keep = Prevalences > Prev_thresh
ibd_ll_f = ibd_ll_data[Vector_keep.index[Vector_keep]]

#Library filter (optional)
ibd_ll_f = ibd_ll_f.filter(regex='agilent', axis=1) #Selection of probes belonging to agilent library

#Include 323 extra probes that are >5% in LL but <5% in IBD
df_extraCC = pd.read_excel("datasetspandas/List_antibodies_to_add_CC.xlsx")

#Concatenate extra peptides with existing dataframe
CCdf = pd.concat([ibd_ll_f, df_extraCC], axis=1)

#Import list of abs to remove and exclude random peptide selection from identical sequences)
df_antibodies_to_remove = pd.read_excel("datasetspandas/Probes_to_remove_CC.xlsx")
abs_list = df_antibodies_to_remove['Name_probe'].values
CCdf.drop(abs_list, axis=1, inplace=True)

#Exclude coagulation-associated peptides
df_serumplasma_antibodies_remove = pd.read_excel("datasetspandas/serumplasmaprobes.xlsx")
sp_list = df_serumplasma_antibodies_remove['Name_probe'].values
CCdf.drop(sp_list, axis=1, inplace=True)

#Rename ID column to allow merging with cohort file
ibd_ll_f = ibd_ll_f.reset_index()
ibd_ll_f.rename(columns={'Unnamed: 0': 'sample_id'}, inplace=True)
ibd_ll_f

#4 Load cohort file
IBDvsHC = pd.read_excel("datasetspandas/IBDvsHC.xlsx")
IBDvsHC

#or: UCvsHC = pd.read_excel("Covariates/UCvsHC-PredictionModel.xlsx") // CDvsHC = pd.read_excel("Covariates/CDvsHC-PredictionModel.xlsx")

#5 Ensure numericity of data input and set indices
ibd_ll_f.apply(pd.to_numeric, errors='ignore')
IBDvsHC.apply(pd.to_numeric, errors='ignore')
IBDvsHC = IBDvsHC.set_index('sample_id')
ibd_ll_f = ibd_ll_f.set_index('sample_id')

#6 Merge datasets to get equal sample numbers
datamerge = IBDvsHC.merge(ibd_ll_f, on=["sample_id"])
datamerge

#7 Fetch cohort classification and remove it from antibody dataframe
df_y = datamerge[['Cohort']]
del datamerge['Cohort']

#8 Import packages for classification analysis and ROC curves
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import auc
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import plot_precision_recall_curve

#Define training and test datasets, 80% train set, 20% test set
trainX, testX, trainy, testy = train_test_split(datamergeUC, df_y, test_size=0.25, random_state=1)

#Fit classification method (RF and logistic regression)
model = LogisticRegression(solver='sag', max_iter=400, penalty='L1/L2', C=) ##determine C-parameter with grid search CV (see below)
model.fit(trainX, trainy.values.ravel())

#Calculate predicted probabilities and predict classes for IBD vs. controls by running model in the test set
probs = model.predict_proba(testX)
preds = probs[:, 1]
y_pred = model.predict(testX)

#Precision and Recall
precision, recall, _ = precision_recall_curve(testy, preds)
f1, auc = f1_score(testy, y_pred), auc(recall, precision)

#summarize scores
print('Logistic: f1=%.3f auc=%.3f' % (f1, auc))

#Precision/recall curve
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set_style("white")
plt.plot(recall, precision, color='blue', label='CD vs HC')
plt.fill_between(recall, precision, color='red', alpha=0.1)
plt.xlim(0, 1.0)
plt.ylim(0, 1.01)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.legend(loc="lower right")
plt.title("Precision - Recall curve (F-score: x.xx, AUC=x.xx)")
plt.savefig("datasetspandas/PR-curve-CDvsHC.jpg", dpi=600, bbox_inches='tight')

#Confusion matrices

#Calculate true negatives (tn), false positives (fp), false negatives (fn) and true positives (tp)
tn, fp, fn, tp = confusion_matrix(testy, y_pred).ravel()
print(f'True Positives: {tp}')
print(f'False Positives: {fp}')
print(f'True Negatives: {tn}')
print(f'False Negatives: {fn}')

#Compute cf matrix
cf_matrix = confusion_matrix(testy, y_pred)

#CF matrix heatmap
x_axis_labels = ['Predicted HC', 'Predicted UC']
y_axis_labels = ['Actual HC', 'Actual UC']
group_names = ['True HC','False UC','False HC','True UC']
group_counts = ["n={0:0.0f}".format(value) for value in cf_matrix.flatten()]
group_percentages = ["{0:.1%}".format(value) for value in cf_matrix.flatten()/np.sum(cf_matrix)]
labels = [f"{v1}\n{v2}\n{v3}" for v1, v2, v3 in zip(group_names,group_counts,group_percentages)]
labels = np.asarray(labels).reshape(2,2)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
heatmap = sns.heatmap(cf_matrix, annot=labels, fmt='', square=True, cmap='Blues', cbar_kws={'label': 'Number of cases in test set (25%)', 'orientation': 'vertical'})
heatmap.set_title("Confusion Matrix: UC vs. HC (test set)", fontdict={'fontsize':12}, fontweight='bold', pad=12)
heatmap.set_yticklabels(labels=y_axis_labels, fontsize=12, rotation=45)
heatmap.set_xticklabels(labels=x_axis_labels, fontsize=12, rotation=45)
heatmap.figure.savefig("datasetspandas/ConfusionMatrix-UCvsHC-logregmodel.jpg", dpi=600, bbox_inches='tight')

#Calculate AUC of predicted probabilities on test set
auc = roc_auc_score(testy, preds)
print('AUC: %.2f' % auc)
print("standard deviation of %0.2f" % (auc.std()))

#Calculate cross-validated AUC
scores = cross_val_score(model, trainX, trainy.values.ravel(), cv=5)
print("%0.2f accuracy with a standard deviation of %0.2f" % (scores.mean(), scores.std()))

#Optimization C parameter using grid search
def Do_gridsearch(trainX,trainy, n_folds):
        #Set values of the grid search; C: float, default=1.0. Inverse of regularization strength, the smaller the C the stronger the regularization. 
        C_values = [0.001, 0.01, 0.1, 1, 10, 100, 1000]
        C_grid = {'C': C_values}
        Logistic_model = LogisticRegression(solver='sag', max_iter=400)
        grid_logReg = GridSearchCV(Logistic_model, C_grid, cv=n_folds, refit=True,n_jobs= -1) #scoring: method to score, n_jobs: jobs in paralel
        grid_logReg.fit(trainX,trainy.values.ravel())

        CV_results = pd.DataFrame( list(zip(C_values, grid_logReg.cv_results_['mean_test_score'])), columns=["C" , "CV"])
        CV_results = CV_results.sort_values(by="CV",ascending=False)

        best_logReg = grid_logReg.best_estimator_
        Accuracy =  best_logReg.score(trainX,trainy)
        return(Accuracy, best_logReg)

Do_gridsearch(trainX,trainy, 10)

#Define fpr, tpr and thresholds for ROC curve
fpr, tpr, thresholds = roc_curve(testy, preds)

#Plotting
def plot_roc_curve(fpr, tpr):
    sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
    sns.set_style("white")
    plt.rcParams["figure.figsize"] = (6,6)
    plt.plot(fpr, tpr, color='orange', label="AUC = 0.91")
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.fill_between(fpr, tpr, color='orange', alpha=0.1)
    plt.title('Antibody repertoire')
    plt.legend(loc="lower right")
    plt.show()

plot_roc_curve(fpr, tpr)
