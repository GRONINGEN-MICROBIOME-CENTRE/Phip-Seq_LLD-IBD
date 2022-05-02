source activate Genomics #conda environment with MEME

Cluster_N=$1

FASTA=../Clusters/$Cluster_N\.fa #FASTA file should be ungapped
Out=MEME_$Cluster_N 


#Ej. script from website
meme $FASTA -protein  -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 3 -maxw 20 -objfun classic -markov_order 0 -oc $Out -seed 99 -minsites 7

