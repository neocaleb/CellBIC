function [exclusivnessMat,cell_grouping]=exclusivnessCal(data_grouping_genecluster,data_grouping_genecluster1,data_grouping_hamming1genecluster_gene,winSizeRatio)
[~,data_grouping_geneclusterSumsortIndex]=sort(sum(data_grouping_genecluster));
[~,sortIndex]=sort(data_grouping_hamming1genecluster_gene(1,:));

gene_group1=data_grouping_hamming1genecluster_gene(sortIndex(1),:)<0.5&~data_grouping_hamming1genecluster_gene(sortIndex(end),:)<0.5;
gene_group2=~data_grouping_hamming1genecluster_gene(sortIndex(1),:)<0.5&data_grouping_hamming1genecluster_gene(sortIndex(end),:)<0.5;
temp2=bwmorph(bwmorph(data_grouping_genecluster1(sortIndex,data_grouping_geneclusterSumsortIndex),'fill'),'clean');
temp3=temp2;
temp3(gene_group2(sortIndex),:)=~temp2(gene_group2(sortIndex),:);

winSize=round(size(data_grouping_genecluster,2)/winSizeRatio);
if winSize==1    
    temp1=smooth(sum(temp3)/size(data_grouping_genecluster,1),winSize+1,'lowess');
    temp1movstd=movstd(temp1,winSize+1);
    [~,peaksLoc]=findpeaks(temp1movstd);
else
    temp1=smooth(sum(temp3)/size(data_grouping_genecluster,1),winSize,'lowess');
    temp1movstd=movstd(temp1,winSize);
    temp1movstd(1:winSize+1)=temp1movstd(winSize+1);
    temp1movstd(size(data_grouping_genecluster,2)-winSize:end)=temp1movstd(size(data_grouping_genecluster,2)-winSize);
    [~,peaksLoc]=findpeaks(temp1movstd);
end

if isempty(peaksLoc)
    exclusivnessMat=[0 0;0 0];
    cell_grouping=ones(1,size(data_grouping_genecluster1,2));
else
    exclusivnessMatTotal=cell(size(peaksLoc));
    exclusivnessMatCriterion05Total=zeros(size(peaksLoc));
    for i=1:size(peaksLoc,1)
        
        exclusivness1(1,1)=mean(mean(temp2(gene_group1(sortIndex),1:peaksLoc(i)-1)));
        exclusivness1(1,2)=1-mean(mean(temp2(gene_group1(sortIndex),peaksLoc(i):end)));
        exclusivness1(2,1)=1-mean(mean(temp2(gene_group2(sortIndex),1:peaksLoc(i)-1)));
        exclusivness1(2,2)=mean(mean(temp2(gene_group2(sortIndex),peaksLoc(i):end)));
        exclusivness2(1,1)=1-mean(mean(temp2(gene_group1(sortIndex),1:peaksLoc(i)-1)));
        exclusivness2(1,2)=mean(mean(temp2(gene_group1(sortIndex),peaksLoc(i):end)));
        exclusivness2(2,1)=mean(mean(temp2(gene_group2(sortIndex),1:peaksLoc(i)-1)));
        exclusivness2(2,2)=1-mean(mean(temp2(gene_group2(sortIndex),peaksLoc(i):end)));
        if sum(sum(exclusivness1))>sum(sum(exclusivness2))
            exclusivnessMat=exclusivness1;
        else
            exclusivnessMat=exclusivness2;
        end
        exclusivnessMatTotal{i}=exclusivnessMat;
        exclusivnessMatCriterion05Total(i)=sum(sum(exclusivnessMat>0.5))==4;
    end
    exclusivnessMatCriterion05Total=find(exclusivnessMatCriterion05Total);
    if size(exclusivnessMatCriterion05Total)
        [~,minIndex]=min(abs(0.5-temp1(peaksLoc(exclusivnessMatCriterion05Total))));
        minIndex=exclusivnessMatCriterion05Total(minIndex);
        cell_grouping=zeros(1,size(data_grouping_genecluster1,2));
        cell_grouping(data_grouping_geneclusterSumsortIndex(1:peaksLoc(minIndex)-1))=1;
        cell_grouping(data_grouping_geneclusterSumsortIndex(peaksLoc(minIndex):end))=2;
        exclusivnessMat=exclusivnessMatTotal{minIndex};
    else
        [~,minIndex]=min(abs(0.5-temp1(peaksLoc)));
        cell_grouping=zeros(1,size(data_grouping_genecluster1,2));
        cell_grouping(data_grouping_geneclusterSumsortIndex(1:peaksLoc(minIndex)-1))=1;
        cell_grouping(data_grouping_geneclusterSumsortIndex(peaksLoc(minIndex):end))=2;
        exclusivnessMat=exclusivnessMatTotal{minIndex};
    end
end