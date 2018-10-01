%% Data imported as RPNegIDS1/RPNEGIDS2 and RPPosIDS1/RPPosIDS2 (Norm/Raw)
RPPosNorm=RPPosIDS1;
RPPosRaw=RPPosIDS2;
clear RPPosIDS1 RPPosIDS2

RPNegNorm=RPNegIDS1;
RPNegRaw=RPNegIDS2;
clear RPNegIDS1 RPNegIDS2

%% Create Structure

RP_POS.DataRaw=RPPosRaw;
RP_POS.Data=RPPosNorm;


%% Remove 'Bad' isotopes dist
j=1;
for i=1:size(RP_POS.Data,1)
    if isnumeric(RP_POS.Data{i,5})
        idx(j)=i;
        j=j+1;
    elseif i>3 && ~ismember(RP_POS.Data{i,5}(1),'1')
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans 
RP_POS.IsotopeIssues.Data=[RP_POS.Data([1:2 6:7],:);RP_POS.Data(idx,:)];
RP_POS.IsotopeIssues.idx=idx;

RP_POS.Data(idx,:)=[];
RP_POS.DataRaw(idx,:)=[];
clear idx

%% Create a grouping list and remove blanks

RP_POS.All_Groups=RP_POS.Data(2,find(ismember(RP_POS.Data(1,:),'Normalised abundance')):end);
for i=2:length(RP_POS.All_Groups)
    if i<7
        RP_POS.All_Groups{i}='Pooled';
    elseif i>=7 && i<12
        RP_POS.All_Groups{i}='Blank';
    elseif i>=12 && i<20
        RP_POS.All_Groups{i}='QC';
    elseif i>=20 && i<96
        RP_POS.All_Groups{i}='Recurr';
    elseif i>=96
        RP_POS.All_Groups{i}='No Recurr';
    end
end

idx=find(ismember(RP_POS.All_Groups,'Blank'));

RP_POS.MASS_RT=cell2mat(RP_POS.Data(4:end,[1 2]))';
RP_POS.All_Titles=RP_POS.Data(3,11:end);
RP_POS.All_Titles(idx)=[];
RP_POS.DataFilt=RP_POS.Data;
RP_POS.DataFilt(:,idx+10)=[];
RP_POS.NormData=cell2mat(RP_POS.Data(4:end,11:end))';
% RP_POS.RawData=cell2mat(RP_POS.DataRaw(4:end,11:end))';
RP_POS.NormData(idx,:)=[];
% RP_POS.RawData(idx,:)=[];

RP_POS.All_Groups(idx)=[];

%% Create Structure for Negative Mode

RP_NEG.DataRaw=RPNegRaw;
RP_NEG.Data=RPNegNorm;


%% Remove 'Bad' isotopes dist
j=1;
for i=1:length(RP_NEG.Data)
    if isnumeric(RP_NEG.Data{i,5})
        idx(j)=i;
        j=j+1;
    elseif i>3 && ~ismember(RP_NEG.Data{i,5}(1),'1')
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans 
RP_NEG.IsotopeIssues.Data=[RP_NEG.Data(1:3,:);RP_NEG.Data(idx,:)];
RP_NEG.IsotopeIssues.idx=idx;

RP_NEG.Data(idx,:)=[];
RP_NEG.DataRaw(idx,:)=[];
clear idx

%% Create a grouping list and remove blanks

RP_NEG.All_Groups=RP_NEG.Data(2,find(ismember(RP_NEG.Data(1,:),'Normalised abundance')):end);
for i=2:length(RP_NEG.All_Groups)
    if i<12
        RP_NEG.All_Groups{i}='Blank';
    elseif i>=12 && i<16
        RP_NEG.All_Groups{i}='Pooled';
    elseif i>=16 && i<24
        RP_NEG.All_Groups{i}='QC';
    elseif i>=24 && i<98
        RP_NEG.All_Groups{i}='Recurr';
    elseif i>=98
        RP_NEG.All_Groups{i}='No Recurr';
    end
end

idx=find(ismember(RP_NEG.All_Groups,'Blank'));

RP_NEG.MASS_RT=cell2mat(RP_NEG.Data(4:end,[1 2]))';
RP_NEG.All_Titles=RP_NEG.Data(3,12:end);
RP_NEG.All_Titles(idx)=[];
RP_NEG.DataFilt=RP_NEG.Data;
RP_NEG.DataFilt(:,idx+11)=[];

RP_NEG.NormData=cell2mat(RP_NEG.Data(4:end,12:end))';
% RP_NEG.RawData=cell2mat(RP_NEG.DataRaw(4:end,12:end))';
RP_NEG.NormData(idx,:)=[];
% RP_NEG.RawData(idx,:)=[];

RP_NEG.All_Groups(idx)=[];

%% Plot Normalized
%Progenesis
normcheckSubplot(RP_NEG.RawData,RP_NEG.NormData,RP_NEG.All_Titles,'one','No Normalization','two','Progenesis Normalization');
normcheckSubplot(RP_POS.RawData,RP_POS.NormData,RP_POS.All_Titles,'one','No Normalization','two','Progenesis Normalization');

%% Boxplot to Visualize

maboxplot(log((RP_NEG.RawData./repmat(median(RP_NEG.RawData),[size(RP_NEG.RawData,1),1]))'),'title','Raw Data','orientation','horizontal')
ylabel('log(Feature/Median Feature Value)')
xlabel('Sample')
set(gca,'XTickLabel',RP_NEG.All_Titles)

maboxplot(log((RP_NEG.NormData./repmat(median(RP_NEG.NormData),[size(RP_NEG.NormData,1),1]))'),'title','Progenesis','orientation','horizontal')
ylabel('log(Feature/Median Feature Value)')
xlabel('Sample')
set(gca,'XTickLabel',RP_NEG.All_Titles)

%% Boxplot to Visualize

maboxplot(log((RP_POS.RawData./repmat(median(RP_POS.RawData),[size(RP_POS.RawData,1),1]))'),'title','Raw Data','orientation','horizontal')
ylabel('log(Feature/Median Feature Value)')
xlabel('Sample')
set(gca,'XTickLabel',RP_POS.All_Titles)

maboxplot(log((RP_POS.MNP./repmat(median(RP_POS.MNP),[size(RP_POS.MNP,1),1]))'),'title','PQN','orientation','horizontal')
ylabel('log(Feature/Median Feature Value)')
xlabel('Sample')
set(gca,'XTickLabel',RP_POS.All_Titles)

%% Average Replicates Negative
%Step 1: Create a matrix with the indeces corresponding to the replicates

AvgInd=zeros(95,2);

for i=1:12 %QC and Pooled
AvgInd(i,1:2)=[i 0];
end
RP_NEG.All_Titles=strrep(RP_NEG.All_Titles, '_1', '');
RP_NEG.All_Titles=strrep(RP_NEG.All_Titles, '_2', '');


j=1; %Samples
for i=13:95
    %disp(j)
    if ~isempty(find(ismember(RP_NEG.All_Titles,['S' num2str(j)])));
        AvgInd(i,:)=find(ismember(RP_NEG.All_Titles,['S' num2str(j)]));
    else
         AvgInd(i,1:2)=0;         
    end
    j=j+1;
end
r=find(sum(AvgInd,2)==0)
AvgInd(r,:)=[];

% Check it
 
 for i=1:length(AvgInd)
     if AvgInd(i,2)~=0;
     A(i,1:2)=RP_NEG.All_Titles(AvgInd(i,:));
     else
     A(i,1)=RP_NEG.All_Titles(AvgInd(i,1));   
     end
 end
 clear A r i j

% Step 2: Average along these rows 
 for i=1:length(AvgInd)
     if AvgInd(i,2)~=0
         RP_NEG.AvgNorm(i,:)=mean(RP_NEG.NormData(AvgInd(i,1:2),:),1); 
         RP_NEG.AvgTitles(i)=RP_NEG.All_Titles(AvgInd(i,1));
     else
         RP_NEG.AvgNorm(i,:)=RP_NEG.NormData(AvgInd(i,1),:); 
         RP_NEG.AvgTitles(i)=RP_NEG.All_Titles(AvgInd(i,1));
     end
 end
 
 RP_NEG.AvgInd=AvgInd;
 RP_NEG.AvgGroups=RP_NEG.All_Groups(RP_NEG.AvgInd(:,1));
 clear A r i j AvgInd
 
 %% Average Replicates Positive
%Step 1: Create a matrix with the indeces corresponding to the replicates
AvgInd=zeros(97,2);

for i=1:14 %QC and Pooled
AvgInd(i,1:2)=[i 0];
end
RP_POS.All_Titles=strrep(RP_POS.All_Titles, '_1', '');
RP_POS.All_Titles=strrep(RP_POS.All_Titles, '_2', '');

j=1; %Samples
for i=15:97
    %disp(j)
    if ~isempty(find(ismember(RP_POS.All_Titles,['S' num2str(j)])));
        AvgInd(i,:)=find(ismember(RP_POS.All_Titles,['S' num2str(j)]));
    else
         AvgInd(i,1:2)=0;         
    end
    j=j+1;
end
[r,~]=find(sum(AvgInd,2)==0);
AvgInd(r,:)=[];
% Check it
 
 for i=1:length(AvgInd)
     if AvgInd(i,2)~=0;
     A(i,1:2)=RP_POS.All_Titles(AvgInd(i,:));
     else
     A(i,1)=RP_POS.All_Titles(AvgInd(i,1));   
     end
 end
 clear A r i j ans

% Step 2: Average along these rows 
 for i=1:length(AvgInd)
     if AvgInd(i,2)~=0
         RP_POS.AvgNorm(i,:)=mean(RP_POS.NormData(AvgInd(i,:),:),1); 
         RP_POS.AvgTitles(i)=RP_POS.All_Titles(AvgInd(i,1));
     else
         RP_POS.AvgNorm(i,:)=RP_POS.NormData(AvgInd(i,1),:); 
         RP_POS.AvgTitles(i)=RP_POS.All_Titles(AvgInd(i,1));
     end
 end
 
 RP_POS.AvgInd=AvgInd;
 RP_POS.AvgGroups=RP_POS.All_Groups(RP_POS.AvgInd(:,1));

 clear A r i j AvgInd idx
 %% Boxplot to Visualize
normcheckSubplot(RP_NEG.NormData,RP_NEG.AvgNorm,RP_NEG.AvgTitles,'one','Normalization','two','Average Normalization');
normcheckSubplot(RP_POS.NormData,RP_POS.AvgNorm,RP_POS.AvgTitles,'one','Normalization','two','Average Normalization');

maboxplot(log((RP_POS.AvgNorm./repmat(median(RP_POS.AvgNorm),[size(RP_POS.AvgNorm,1),1]))'),'title','Averaged Progenesis Normalized Data','orientation','horizontal')
ylabel('log(Feature/Median Feature Value)')
xlabel('Sample')
set(gca,'XTickLabel',RP_POS.AvgTitles)

maboxplot(log((RP_NEG.AvgNorm./repmat(median(RP_NEG.AvgNorm),[size(RP_NEG.AvgNorm,1),1]))'),'title','Averaged Progenesis Normalized Data','orientation','horizontal')
ylabel('log(Feature/Median Feature Value)')
xlabel('Sample')
set(gca,'XTickLabel',RP_NEG.AvgTitles)


%% Remove Pooled, QC, Sample 77 
RP_NEG.AvgNorm_wo77=RP_NEG.AvgNorm([13:85 87:end],:);
RP_NEG.AvgTitles_wo77=RP_NEG.AvgTitles([13:85 87:end]);
RP_NEG.AvgGroups_wo77=RP_NEG.AvgGroups([13:85 87:end]);

%% Remove Pooled, QC, Sample 77 
RP_POS.AvgNorm_wo77=RP_POS.AvgNorm([15:87 89:end],:);
RP_POS.AvgTitles_wo77=RP_POS.AvgTitles([15:87 89:end]);
RP_POS.AvgGroups_wo77=RP_POS.AvgGroups([15:87 89:end]);

%% Create Y Vector Positive
RP_POS.Y=zeros(length(RP_POS.AvgGroups_wo77),1);
RP_POS.Y(find(ismember(RP_POS.AvgGroups_wo77,'No Recurr')))=2;  %2 is No Reccurance, 1 is Recurrance
RP_POS.Y(find(ismember(RP_POS.AvgGroups_wo77,'Recurr')))=1;

RP_POS.ID=RP_POS.DataFilt(3:end,[6 7 9 10 4]);
%% Get the significance of each Data set with an FDR correction
%RP POS
norecurr=RP_POS.AvgNorm_wo77(RP_POS.Y==2,:); %no recur
recurr=RP_POS.AvgNorm_wo77(RP_POS.Y==1,:); %recur
norecurr_tit=RP_POS.AvgTitles_wo77(RP_POS.Y==2); %no recur
recurr_tit=RP_POS.AvgTitles_wo77(RP_POS.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

for i=1:length(norecurr)
    [~,RP_POS.Sign.Pvalues(i)]=ttest2(norecurr(:,i),recurr(:,i),'Vartype','unequal');
end

[~, RP_POS.AdjSign.Crit_pvalue, RP_POS.AdjSign.Adjpvalues]=fdr_bh(RP_POS.Sign.Pvalues); %if using FDR 
RP_POS.Sign.Signindx=find(RP_POS.Sign.Pvalues<0.05);
RP_POS.Sign.SignInt=RP_POS.AvgNorm_wo77(:,RP_POS.Sign.Signindx);
RP_POS.Sign.SignMRT=RP_POS.MASS_RT(:,RP_POS.Sign.Signindx)';
RP_POS.Sign.SignID=RP_POS.ID([1 RP_POS.Sign.Signindx],:);

RP_POS.AdjSign.AdjSignindx=find(RP_POS.AdjSign.Adjpvalues<0.05);
RP_POS.AdjSign.AdjSignInt=RP_POS.AvgNorm_wo77(:,RP_POS.AdjSign.AdjSignindx);
RP_POS.AdjSign.AdjSignMRT=RP_POS.MASS_RT(:,RP_POS.AdjSign.AdjSignindx)';
RP_POS.AdjSign.AdjSignID=RP_POS.ID([1 RP_POS.AdjSign.AdjSignindx+1],:);


RP_POS.Sign.SignFoldChange=VolcanoPlot(norecurr,recurr,RP_POS.Sign.Pvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',RP_POS.MASS_RT(1,:),'Features2',RP_POS.MASS_RT(2,:),'Groups',Groups);
RP_POS.Sign.SignFoldChange.ID=RP_POS.ID([1 RP_POS.Sign.SignFoldChange.idx+1],:);

RP_POS.AdjSign.AdjFoldChange=VolcanoPlot(norecurr,recurr,RP_POS.AdjSign.Adjpvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',RP_POS.MASS_RT(1,:),'Features2',RP_POS.MASS_RT(2,:),'Groups',Groups);
RP_POS.AdjSign.AdjFoldChange.ID=RP_POS.ID([1 RP_POS.AdjSign.AdjFoldChange.idx+1],:);

clear norecurr reccur norecurr_tit reccur_tit Groups i

%% Create Y Vector Negative
RP_NEG.Y=zeros(length(RP_NEG.AvgGroups_wo77),1);
RP_NEG.Y(find(ismember(RP_NEG.AvgGroups_wo77,'No Recurr')))=2;  %2 is No Reccurance, 1 is Recurrance
RP_NEG.Y(find(ismember(RP_NEG.AvgGroups_wo77,'Recurr')))=1;

RP_NEG.ID=RP_NEG.DataFilt(3:end,[6 7 9 10 4]);
%% Get the significance of each Data set with an FDR correction
%RP NEG
norecurr=RP_NEG.AvgNorm_wo77(RP_NEG.Y==2,:); %no recur
recurr=RP_NEG.AvgNorm_wo77(RP_NEG.Y==1,:); %recur
norecurr_tit=RP_NEG.AvgTitles_wo77(RP_NEG.Y==2); %no recur
recurr_tit=RP_NEG.AvgTitles_wo77(RP_NEG.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

for i=1:length(norecurr)
    [~,RP_NEG.Sign.Pvalues(i)]=ttest2(norecurr(:,i),recurr(:,i),'Vartype','unequal');
end

[~, RP_NEG.AdjSign.Crit_pvalue, RP_NEG.AdjSign.Adjpvalues]=fdr_bh(RP_NEG.Sign.Pvalues); %if using FDR 
RP_NEG.Sign.Signindx=find(RP_NEG.Sign.Pvalues<0.05);
RP_NEG.Sign.SignInt=RP_NEG.AvgNorm_wo77(:,RP_NEG.Sign.Signindx);
RP_NEG.Sign.SignMRT=RP_NEG.MASS_RT(:,RP_NEG.Sign.Signindx)';
RP_NEG.Sign.SignID=RP_NEG.ID([1 RP_NEG.Sign.Signindx+1],:);

RP_NEG.AdjSign.AdjSignindx=find(RP_NEG.AdjSign.Adjpvalues<0.05);
RP_NEG.AdjSign.AdjSignInt=RP_NEG.AvgNorm_wo77(:,RP_NEG.AdjSign.AdjSignindx);
RP_NEG.AdjSign.AdjSignMRT=RP_NEG.MASS_RT(:,RP_NEG.AdjSign.AdjSignindx)';
RP_NEG.AdjSign.AdjSignID=RP_NEG.ID([1 RP_NEG.AdjSign.AdjSignindx+1],:);


RP_NEG.Sign.SignFoldChange=VolcanoPlot(norecurr,recurr,RP_NEG.Sign.Pvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',RP_NEG.MASS_RT(1,:),'Features2',RP_NEG.MASS_RT(2,:),'Groups',Groups);
RP_NEG.Sign.SignFoldChange.ID=RP_NEG.ID([1 RP_NEG.Sign.SignFoldChange.idx],:);

RP_NEG.AdjSign.AdjFoldChange=VolcanoPlot(norecurr,recurr,RP_NEG.AdjSign.Adjpvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',RP_NEG.MASS_RT(1,:),'Features2',RP_NEG.MASS_RT(2,:),'Groups',Groups);
RP_NEG.AdjSign.AdjFoldChange.ID=RP_NEG.ID([1 RP_NEG.AdjSign.AdjFoldChange.idx],:);

clear norecurr reccur norecurr_tit reccur_tit Groups i