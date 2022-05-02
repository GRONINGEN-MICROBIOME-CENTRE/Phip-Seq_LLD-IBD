source activate ../../conda_env/ATLAS/

MSA=$1 #Sequences_cluster_smoking.txt
OUT=$2 #Smoking_probes.mat

#MSA=Sequences_cluster.txt
#OUT=Cluster3.mat

clustalo -i $MSA --distmat-out=$OUT --full


