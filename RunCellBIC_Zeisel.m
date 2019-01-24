%%%%%%%%%%%%% preparation of data %%%%%%%%%%%%%
load 'Zeisel.mat'
%%%% Parameter Setting %%%%
MaxCellInCluster=100;
clusterScoreWeight=0.7;
winSizeRatio=15;
minGeneGroupNum=10;
maxNumClust=5;
distanceFromSeedCutoff=0.2;
seedDistanceCutoff=0.7;
%%%% Run CellBIC step1 %%%%
[data_grouping1_iter,bimodal_gene_list_iter,genecluster_total_iter,cell_grouping_total_iter,exclusivness_total_iter,log_data_select_iter,iter_depth]=CellBIC_step1(log_data,MaxCellInCluster,clusterScoreWeight,winSizeRatio,minGeneGroupNum,maxNumClust,distanceFromSeedCutoff,seedDistanceCutoff);
%%%% Run CellBIC step2 %%%%
numClust=7;
[cell_grouping_alter,clusterScore_alter]=CellBIC_step2(genecluster_total_iter,cell_grouping_total_iter,exclusivness_total_iter,log_data_select_iter,iter_depth,numClust);
