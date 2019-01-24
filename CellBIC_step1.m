function [data_grouping1_iter,bimodal_gene_list_iter,genecluster_total_iter,cell_grouping_total_iter,exclusivness_total_iter,log_data_select_iter,iter_depth]=CellBIC_step1(log_data,MaxCellInCluster,clusterScoreWeight,winSizeRatio,minGeneGroupNum,maxNumClust,distanceFromSeedCutoff,seedDistanceCutoff)
iter=1;
[data_grouping1,bimodal_gene_list,genecluster_total,cell_grouping_total,exclusivness_total]=bimodal_clustering(log_data,winSizeRatio,minGeneGroupNum,maxNumClust,distanceFromSeedCutoff,seedDistanceCutoff);
data_grouping1_iter{1}=data_grouping1;
bimodal_gene_list_iter{1}=bimodal_gene_list;
genecluster_total_iter{1}=genecluster_total;
cell_grouping_total_iter{1}=cell_grouping_total;
exclusivness_total_iter{1}=exclusivness_total;
log_data_select_iter{1}=find(ones(1,size(log_data,2)));

iter_did=[1];
iter_mother=[0];
iter_mother_clustering=[0];
iter_mother_cluster=[0];
iter_depth=[1];
[~,genecluster_selectIndex]=max((exclusivness_total-clusterScoreWeight).*sum(genecluster_total));
cell_grouping=cell_grouping_total(genecluster_selectIndex,:);
for i=1:max(cell_grouping)
    if sum(cell_grouping==i)>MaxCellInCluster
        iter=iter+1;
        iter_did(iter)=0;
        iter_mother(iter)=1;
        iter_mother_clustering(iter)=genecluster_selectIndex;
        iter_mother_cluster(iter)=i;
        iter_depth(iter)=iter_depth(iter_mother(iter))+1;
        log_data_select_iter{iter}=log_data_select_iter{iter_mother(iter)}(cell_grouping==i);        
    end
end

while sum(iter_did==0)>0
    for iter=1:size(iter_did,2)
        if iter_did(iter)==0            
            [data_grouping1,bimodal_gene_list,genecluster_total,cell_grouping_total,exclusivness_total]=bimodal_clustering(log_data(:,log_data_select_iter{iter}),winSizeRatio,minGeneGroupNum,maxNumClust,distanceFromSeedCutoff,seedDistanceCutoff);
            data_grouping1_iter{iter}=data_grouping1;
            bimodal_gene_list_iter{iter}=bimodal_gene_list;            
            genecluster_total_iter{iter}=genecluster_total;
            cell_grouping_total_iter{iter}=cell_grouping_total;
            exclusivness_total_iter{iter}=exclusivness_total;
            iter_did(iter)=1;
            
            [~,genecluster_selectIndex]=max((exclusivness_total-clusterScoreWeight).*sum(genecluster_total));
            cell_grouping=cell_grouping_total(genecluster_selectIndex,:);
            iter_previous=iter;
            iter=size(iter_did,2);
            for i=1:max(cell_grouping)
                if sum(cell_grouping==i)>MaxCellInCluster
                    iter=iter+1;
                    iter_did(iter)=0;
                    iter_mother(iter)=iter_previous;
                    iter_mother_clustering(iter)=genecluster_selectIndex;
                    iter_mother_cluster(iter)=i;
                    iter_depth(iter)=iter_depth(iter_mother(iter))+1;
                    log_data_select_iter{iter}=log_data_select_iter{iter_mother(iter)}(cell_grouping==i);
                end
            end
        end
    end
end