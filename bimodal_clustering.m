function [data_grouping1,bimodal_gene_list,genecluster_total,cell_grouping_total,exclusivness_total]=bimodal_clustering(log_data_temp,winSizeRatio,minGeneGroupNum,maxNumClust,distanceFromSeedCutoff,seedDistanceCutoff)
%%%%%%%%%%%%% calculation of biomodal metric %%%%%%%%%%%%%
geneFilterIndex=find(quantile(log_data_temp,0.95,2)>1&std(log_data_temp')'>0.001);
options = statset('MaxIter',10000);
bimodalMetric=zeros(size(geneFilterIndex,1),6);
noise_number=10;
randn('seed',1);
for geneIndex = 1:size(geneFilterIndex,1)
    data_temp1=log_data_temp(geneFilterIndex(geneIndex),:)';    
    data_temp1=repmat(data_temp1,noise_number,1)+abs(randn(size(data_temp1,1)*noise_number,1))*std(data_temp1)/sqrt(size(data_temp1,1))/10;
    
    obj1=gmdistribution.fit(data_temp1,2,'Options',options);
    if obj1.mu(1)<=obj1.mu(2)
        bimodalMetric(geneIndex,1:2)=obj1.mu';
        bimodalMetric(geneIndex,3)=sqrt(obj1.Sigma(1));
        bimodalMetric(geneIndex,4)=sqrt(obj1.Sigma(2));
        bimodalMetric(geneIndex,5:6)=obj1.PComponents;
    else
        bimodalMetric(geneIndex,2)=obj1.mu(1);
        bimodalMetric(geneIndex,1)=obj1.mu(2);
        bimodalMetric(geneIndex,4)=sqrt(obj1.Sigma(1));
        bimodalMetric(geneIndex,3)=sqrt(obj1.Sigma(2));
        bimodalMetric(geneIndex,6)=obj1.PComponents(1);
        bimodalMetric(geneIndex,5)=obj1.PComponents(2);
    end
end
data_grouping1=zeros(size(log_data_temp(geneFilterIndex,:)));
data_bimodal_tstatistics=zeros(size(geneFilterIndex,1),1);
for geneIndex = 1:size(geneFilterIndex,1)
    data_temp_zscore1=abs(log_data_temp(geneFilterIndex(geneIndex),:)-bimodalMetric(geneIndex,1))/bimodalMetric(geneIndex,3);
    data_temp_zscore2=abs(log_data_temp(geneFilterIndex(geneIndex),:)-bimodalMetric(geneIndex,2))/bimodalMetric(geneIndex,4);
    data_grouping1(geneIndex,:)=(data_temp_zscore1-data_temp_zscore2)>0;
    data_grouping1(geneIndex,log_data_temp(geneFilterIndex(geneIndex),:)<bimodalMetric(geneIndex,1))=0;
    data_grouping1(geneIndex,log_data_temp(geneFilterIndex(geneIndex),:)>bimodalMetric(geneIndex,2))=1;
    data_bimodal_se=sqrt((bimodalMetric(geneIndex,4)^2/bimodalMetric(geneIndex,6)+bimodalMetric(geneIndex,3)^2/bimodalMetric(geneIndex,5))/size(log_data_temp,2));
    data_bimodal_tstatistics(geneIndex)=(bimodalMetric(geneIndex,2)-bimodalMetric(geneIndex,1))/data_bimodal_se;
end

%%%%%%%%%%%%% calculating hamming distance %%%%%%%%%%%%%
bimodalCriterion=bimodalMetric(:,5)>0.1&bimodalMetric(:,6)>0.1&data_bimodal_tstatistics>10&bimodalMetric(:,2)>1;
bimodal_gene_list=geneFilterIndex(bimodalCriterion);
data_grouping1=data_grouping1(bimodalCriterion,:);
data_grouping_hamming1=pdist(data_grouping1,'hamming');
data_grouping_hamming=data_grouping_hamming1;
data_grouping_hamming(data_grouping_hamming>0.5)=1-data_grouping_hamming(data_grouping_hamming>0.5);
data_grouping_hamming_square1=squareform(data_grouping_hamming1);
data_grouping_hamming_square=squareform(data_grouping_hamming);
%%%%%%%%%%%%% new clustering %%%%%%%%%%%%%
[data_grouping_hamming1sort,data_grouping_hamming1sortIndex]=sort(data_grouping_hamming1,'Descend');
IndexTemp=zeros(size(data_grouping_hamming1,2),2);
IndexCount=1;
for i=1:size(bimodal_gene_list,1)-1
    for j=2:size(bimodal_gene_list,1)
        if i<j
            IndexTemp(IndexCount,:)=[j,i];
            IndexCount=IndexCount+1;
        end
    end
end
if sum(data_grouping_hamming1sort>seedDistanceCutoff)/size(data_grouping_hamming1sort,2)>0.001
    data_grouping_hamming_square1sortIndex=IndexTemp(data_grouping_hamming1sortIndex(data_grouping_hamming1sort>seedDistanceCutoff),:);
else
    data_grouping_hamming_square1sortIndex=IndexTemp(data_grouping_hamming1sortIndex(data_grouping_hamming1sort>0.6),:);
end

genecluster_total=[];
exclusivness_total=[];
cell_grouping_total=[];
gene_grouping_total=[];
genecluster_check=zeros(size(bimodal_gene_list,1),1);
seedIndex=1;
for seedIndex=1:size(data_grouping_hamming_square1sortIndex,1)
    if sum(data_grouping_hamming_square1sortIndex(seedIndex,1)==find(genecluster_check))+sum(data_grouping_hamming_square1sortIndex(seedIndex,2)==find(genecluster_check))==0
        genecluster=data_grouping_hamming_square1sortIndex(seedIndex,:);
        distanceFromSeed=mean(mean(data_grouping_hamming_square(data_grouping_hamming_square1sortIndex(seedIndex,:),genecluster)));        
        data_grouping_genecluster=data_grouping1(genecluster(1),:);
        data_grouping_genecluster=[data_grouping_genecluster;~data_grouping1(genecluster(2),:)];        
        [~,data_grouping_hamming_seed_sortIndex]=sort(mean(data_grouping_hamming_square(data_grouping_hamming_square1sortIndex(seedIndex,:),:)));        
        
        for j=1:size(data_grouping_hamming_seed_sortIndex,2)
            if ~sum(data_grouping_hamming_seed_sortIndex(j)==genecluster) && ~sum(data_grouping_hamming_seed_sortIndex(j)==find(genecluster_check))                
                genecluster=[genecluster,data_grouping_hamming_seed_sortIndex(j)];
                distanceFromSeed=[distanceFromSeed,mean(mean(data_grouping_hamming_square(data_grouping_hamming_square1sortIndex(seedIndex,:),genecluster)))];
                               
                
                if data_grouping_hamming_square1(genecluster(1),genecluster(end))<0.5
                    data_grouping_genecluster=[data_grouping_genecluster;data_grouping1(genecluster(end),:)];
                else
                    data_grouping_genecluster=[data_grouping_genecluster;~data_grouping1(genecluster(end),:)];
                end                
                
                if distanceFromSeed(end)>=distanceFromSeedCutoff
                    break
                end
            end
        end
        genecluster_check(genecluster)=1;
        if size(genecluster,2)>=minGeneGroupNum
            data_grouping_genecluster1=data_grouping1(genecluster,:);
            data_grouping_hamming1genecluster_gene=data_grouping_hamming_square1(genecluster,genecluster);
            [exclusivnessMat,cell_grouping]=exclusivnessCal(data_grouping_genecluster,data_grouping_genecluster1,data_grouping_hamming1genecluster_gene,winSizeRatio);
                        
            if sum(sum(exclusivnessMat>0.5))==4
                genecluster_total=[genecluster_total,zeros(size(bimodal_gene_list))];
                genecluster_total(genecluster,end)=1;
                exclusivness_total=[exclusivness_total,mean(mean(exclusivnessMat))];
                cell_grouping_total=[cell_grouping_total;cell_grouping];
            else
                genecluster_temp=genecluster;
                data_grouping_genecluster_temp=data_grouping_genecluster;
                data_grouping_genecluster_temp1=data_grouping_genecluster1;
                [exclusivnessMatPrev,cell_groupingPrev]=exclusivnessCal(data_grouping_genecluster(1:minGeneGroupNum,:),data_grouping_genecluster1(1:minGeneGroupNum,:),data_grouping_hamming1genecluster_gene(1:minGeneGroupNum,1:minGeneGroupNum),winSizeRatio);
                if sum(sum(exclusivnessMatPrev>0.5))==4
                    for i=minGeneGroupNum+1:size(genecluster,2)
                        [exclusivnessMat,cell_grouping]=exclusivnessCal(data_grouping_genecluster(1:i,:),data_grouping_genecluster1(1:i,:),data_grouping_hamming1genecluster_gene(1:i,1:i),winSizeRatio);
                        if sum(sum(exclusivnessMat>0.5))==4
                            exclusivnessMatPrev=exclusivnessMat;
                            cell_groupingPrev=cell_grouping;
                        else
                            exclusivnessMat=exclusivnessMatPrev;
                            cell_grouping=cell_groupingPrev;
                            genecluster_total=[genecluster_total,zeros(size(bimodal_gene_list))];
                            genecluster_total(genecluster(1:i-1),end)=1;
                            exclusivness_total=[exclusivness_total,mean(mean(exclusivnessMat))];
                            cell_grouping_total=[cell_grouping_total;cell_grouping];
                            break
                        end
                    end
                end
            end
        end
    end
    if size(genecluster_total,2)==maxNumClust
        break
    end
end