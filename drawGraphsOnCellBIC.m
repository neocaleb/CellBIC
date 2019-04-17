function drawGraphsOnCellBIC(data_grouping1_iter,bimodal_gene_list_iter,genecluster_total_iter,cell_grouping_total_iter,log_data_select_iter,iter_depth,cell_id,gene_name,log_data,outputFolder)
for iter=1:size(genecluster_total_iter,2)
    data_grouping1=data_grouping1_iter{iter};
    bimodal_gene_list=bimodal_gene_list_iter{iter};
    genecluster_total=genecluster_total_iter{iter};
    cell_grouping_total=cell_grouping_total_iter{iter};
    log_data_select=log_data_select_iter{iter};
    
    data_grouping_hamming1=pdist(data_grouping1,'hamming');
    data_grouping_hamming_square1=squareform(data_grouping_hamming1);
    
    for clusterIndex=1:size(genecluster_total,2)
        genecluster=find(genecluster_total(:,clusterIndex));
        cell_grouping=cell_grouping_total(clusterIndex,:);
        
        [~,sortIndex]=sort(data_grouping_hamming_square1(genecluster(1),genecluster));
        data_grouping_genecluster=data_grouping1(genecluster(sortIndex(1)),:);
        for i=2:size(sortIndex,2)
            if data_grouping_hamming_square1(genecluster(sortIndex(1)),genecluster(sortIndex(i)))<0.5
                data_grouping_genecluster=[data_grouping_genecluster;data_grouping1(genecluster(sortIndex(i)),:)];
            else
                data_grouping_genecluster=[data_grouping_genecluster;~data_grouping1(genecluster(sortIndex(i)),:)];
            end
        end
        [~,data_grouping_geneclusterSumsortIndex]=sort(sum(data_grouping_genecluster));
        
        data_grouping_hamming1genecluster_gene=data_grouping_hamming_square1(genecluster,genecluster);
        gene_group1=data_grouping_hamming1genecluster_gene(sortIndex(1),:)<0.5&~data_grouping_hamming1genecluster_gene(sortIndex(end),:)<0.5;
        gene_group2=~data_grouping_hamming1genecluster_gene(sortIndex(1),:)<0.5&data_grouping_hamming1genecluster_gene(sortIndex(end),:)<0.5;
        
        f = figure('Visible','off');
        imagesc(data_grouping1(genecluster(sortIndex),data_grouping_geneclusterSumsortIndex))
        yticks([1:size(genecluster,1)])
        yticklabels(gene_name(bimodal_gene_list(genecluster(sortIndex))))
        saveas(f,[outputFolder,'/heatmap.grouping.iter',num2str(iter),'.genecluster',num2str(clusterIndex),'.pdf'])
        
        f = figure('Visible','off');
        imagesc(log_data(bimodal_gene_list(genecluster(sortIndex)),log_data_select(data_grouping_geneclusterSumsortIndex)))
        yticks([1:size(genecluster,1)])
        yticklabels(gene_name(bimodal_gene_list(genecluster(sortIndex))))
        saveas(f,[outputFolder,'/heatmap.logdata.iter',num2str(iter),'.genecluster',num2str(clusterIndex),'.pdf'])
                
        f = figure('Visible','off');
        imagesc(cell_grouping(data_grouping_geneclusterSumsortIndex))
        set(gcf, 'Position', [0 0 560 100]);
        caxis([0 2])
        colormap jet
        saveas(f,[outputFolder,'/Index.cell_grouping.iter',num2str(iter),'.genecluster',num2str(clusterIndex),'.pdf'])
        
        ofile = fopen([outputFolder,'/genelist.iter',num2str(iter),'.genecluster',num2str(clusterIndex),'.txt'],'w');
        for i=1:size(genecluster,1)
            fprintf(ofile,'%s\n',gene_name(bimodal_gene_list(genecluster(i))));
        end
        fclose(ofile);
        
        ofile = fopen([outputFolder,'/cell_id.clustering.iter',num2str(iter),'.genecluster',num2str(clusterIndex),'.txt'],'w');
        for i=1:size(log_data_select,2)
            fprintf(ofile,'%s\t%d\n',cell_id(log_data_select(i)),cell_grouping(i));
        end
        fclose(ofile);
        
        ofile = fopen([outputFolder,'/genelist_group1.iter',num2str(iter),'.genecluster',num2str(clusterIndex),'.txt'],'w');
        for i=1:size(genecluster,1)
            if gene_group1(i)
                fprintf(ofile,'%s\n',gene_name(bimodal_gene_list(genecluster(i))));
            end
        end
        fclose(ofile);
        
        ofile = fopen([outputFolder,'/genelist_group2.iter',num2str(iter),'.genecluster',num2str(clusterIndex),'.txt'],'w');
        for i=1:size(genecluster,1)
            if gene_group2(i)
                fprintf(ofile,'%s\n',gene_name(bimodal_gene_list(genecluster(i))));
            end
        end
        fclose(ofile);
    end
end