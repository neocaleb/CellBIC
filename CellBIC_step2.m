function [cell_grouping_alter,clusterScore_alter]=CellBIC_step2(genecluster_total_iter,cell_grouping_total_iter,exclusivness_total_iter,log_data_select_iter,iter_depth,numClust)
iter=1;
genecluster_total=genecluster_total_iter{iter};
cell_grouping_total=cell_grouping_total_iter{iter};
exclusivness_total=exclusivness_total_iter{iter};
log_data_select=log_data_select_iter{iter};
cell_grouping_alter=cell_grouping_total;
clusterScore_alter=(exclusivness_total-0.7).*log2(sum(genecluster_total));
depth=2;
numClustCurrent=max(max(cell_grouping_alter));
while numClustCurrent<numClust && depth<=max(iter_depth)
    same_depth=0;
    for i=find(iter_depth==depth)
        if ~isempty(cell_grouping_total_iter{i})
            same_depth=same_depth+1;
        end
    end
    for i=1:min(numClust-numClustCurrent,same_depth)
        cell_grouping_alter_previous=cell_grouping_alter;
        clusterScore_alter_previous=clusterScore_alter;
        cell_grouping_alter=[];
        clusterScore_alter=[];
        for alterIndex=1:size(cell_grouping_alter_previous,1)
            for iter=2:size(iter_depth,2)
                if iter_depth(iter)==depth && ~isempty(cell_grouping_total_iter{iter})
                    genecluster_total=genecluster_total_iter{iter};
                    cell_grouping_total=cell_grouping_total_iter{iter};
                    exclusivness_total=exclusivness_total_iter{iter};
                    log_data_select=log_data_select_iter{iter};
                    if size(unique(cell_grouping_alter_previous(alterIndex,log_data_select)),2)==1
                        for clusterIndex=1:size(cell_grouping_total,1)
                            cell_grouping=cell_grouping_alter_previous(alterIndex,:);
                            cell_grouping(cell_grouping>max(cell_grouping(log_data_select)))=cell_grouping(cell_grouping>max(cell_grouping(log_data_select)))+1;
                            cell_grouping(log_data_select(cell_grouping_total(clusterIndex,:)==2))=max(cell_grouping(log_data_select))+1;
                            if isempty(cell_grouping_alter) || sum(sum(cell_grouping_alter==cell_grouping,2)==size(cell_grouping,2))==0
                                cell_grouping_alter=[cell_grouping_alter;cell_grouping];
                                clusterScore_total=(exclusivness_total-0.7).*log2(sum(genecluster_total));
                                clusterScore_alter=[clusterScore_alter (clusterScore_alter_previous(alterIndex)*(numClustCurrent-1)+clusterScore_total(clusterIndex))/numClustCurrent];
                            end
                        end
                    end
                end
            end
        end
        numClustCurrent=max(max(cell_grouping_alter));
    end
    if numClustCurrent<numClust
        depth=depth+1;
    end
end
if numClust>max(max(cell_grouping_alter))
    disp(['Cells cannot be clustered into ',num2str(numClust),' clusters'])
end
