%% Data imported as HILICNegIDS1/HILICNEGIDS2 and HILICPosIDS1/HILICPosIDS2 (Norm/Raw)
HILICPosNorm=HILICPosID2;
HILICPosRaw=HILICPosIDS2;
clear HILICPosIDS1 HILICPosIDS2

HILICNegNorm=HILICNegIDS1;
HILICNegRaw=HILICNegIDS2;
clear HILICNegIDS1 HILICNegIDS2

%% Create Structure

HILIC_POS.DataRaw=HILICPosRaw;
HILIC_POS.Data=HILICPosNorm;

%% Remove 'Bad' isotopes dist
j=1;
for i=1:size(HILIC_POS.Data,1)
    if isnumeric(HILIC_POS.Data{i,5})
        idx(j)=i;
        j=j+1;
    elseif i>3 && ~ismember(HILIC_POS.Data{i,5}(1),'1')
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans 
HILIC_POS.IsotopeIssues.Data=[HILIC_POS.Data(1:3,:);HILIC_POS.Data(idx,:)];
HILIC_POS.IsotopeIssues.idx=idx;

HILIC_POS.Data(idx,:)=[];
% HILIC_POS.DataRaw(idx,:)=[];
clear idx

%% Create a grouping list and remove blanks

HILIC_POS.All_Groups=HILIC_POS.Data(2,find(ismember(HILIC_POS.Data(1,:),'Normalised abundance')):end);
for i=2:length(HILIC_POS.All_Groups)
    if i<26
        HILIC_POS.All_Groups{i}='Blank';
    elseif i>=26 && i<37
        HILIC_POS.All_Groups{i}='Pooled';
    elseif i>=37 && i<39
        HILIC_POS.All_Groups{i}='DDAx';
    elseif i>=39 && i<153
        HILIC_POS.All_Groups{i}='Recurr';
    elseif i>=153
        HILIC_POS.All_Groups{i}='No Recurr';
    end
end
idx=find(ismember(HILIC_POS.All_Groups,'Blank') | ismember(HILIC_POS.All_Groups,'DDAx'));
HILIC_POS.MASS_RT=cell2mat(HILIC_POS.Data(4:end,[1 2]))';
HILIC_POS.All_Titles=HILIC_POS.Data(3,12:end);


HILIC_POS.All_Titles(idx)=[];
HILIC_POS.DataFilt=HILIC_POS.Data;
HILIC_POS.DataFilt(:,idx+10)=[];
HILIC_POS.NormData=cell2mat(HILIC_POS.Data(4:end,12:end))';
% HILIC_POS.RawData=cell2mat(HILIC_POS.DataRaw(4:end,11:end))';
HILIC_POS.NormData(idx,:)=[];
% HILIC_POS.RawData(idx,:)=[];

HILIC_POS.All_Groups(idx)=[];

%% Create Structure for Negative Mode

HILIC_NEG.DataRaw=HILICNegRaw;
HILIC_NEG.Data=HILICNegNorm;


%% Remove 'Bad' isotopes dist
j=1;
for i=1:size(HILIC_NEG.Data,1)
    if isnumeric(HILIC_NEG.Data{i,5})
        idx(j)=i;
        j=j+1;
    elseif i>3 && ~ismember(HILIC_NEG.Data{i,5}(1),'1')
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans 
HILIC_NEG.IsotopeIssues.Data=[HILIC_NEG.Data(1:3,:);HILIC_NEG.Data(idx,:)];
HILIC_NEG.IsotopeIssues.idx=idx;

HILIC_NEG.Data(idx,:)=[];
% HILIC_NEG.DataRaw(idx,:)=[];
clear idx

%% Create a grouping list and remove blanks

HILIC_NEG.All_Groups=HILIC_NEG.Data(2,find(ismember(HILIC_NEG.Data(1,:),'Normalised abundance')):end);
for i=2:length(HILIC_NEG.All_Groups)
    if i<24
        HILIC_NEG.All_Groups{i}='Blank';
    elseif i>=24 && i<36
        HILIC_NEG.All_Groups{i}='Pooled';
    elseif i>=36 && i<59
        HILIC_NEG.All_Groups{i}='QC';
    elseif i>=59 && i<61
        HILIC_NEG.All_Groups{i}='DDAx';
    elseif i>=61 && i<195
        HILIC_NEG.All_Groups{i}='No Recurr';
    elseif i>=195
        HILIC_NEG.All_Groups{i}='Recurr';
    end
end
idx=find(ismember(HILIC_NEG.All_Groups,'Blank') | ismember(HILIC_NEG.All_Groups,'DDAx'));

HILIC_NEG.MASS_RT=cell2mat(HILIC_NEG.Data(4:end,[1 2]))';
HILIC_NEG.All_Titles=HILIC_NEG.Data(3,11:end);
HILIC_NEG.All_Titles(idx)=[];
HILIC_NEG.DataFilt=HILIC_NEG.Data;
HILIC_NEG.DataFilt(:,idx+10)=[];
HILIC_NEG.NormData=cell2mat(HILIC_NEG.Data(4:end,11:end))';
% HILIC_NEG.RawData=cell2mat(HILIC_NEG.DataRaw(4:end,11:end))';
HILIC_NEG.NormData(idx,:)=[];
% HILIC_NEG.RawData(idx,:)=[];

HILIC_NEG.All_Groups(idx)=[];

clear i idx
%% Plot Normalized
%Progenesis
normcheckSubplot(HILIC_NEG.RawData,HILIC_NEG.NormData,HILIC_NEG.All_Titles,'one','No Normalization','two','Progenesis Normalization');
normcheckSubplot(HILIC_POS.RawData,HILIC_POS.NormData,HILIC_POS.All_Titles,'one','No Normalization','two','Progenesis Normalization');

%% Average Replicates Negative
%Step 1: Create a matrix with the indeces corresponding to the replicates

AvgInd=zeros(121,3);

for i=1:34 % Pooled
AvgInd(i,1:3)=[i 0 0];
end

HILIC_NEG.All_Titles=strrep(HILIC_NEG.All_Titles, '_1', '');
HILIC_NEG.All_Titles=strrep(HILIC_NEG.All_Titles, '_2', '');
HILIC_NEG.All_Titles=strrep(HILIC_NEG.All_Titles, '_3', '');
%%
j=1; %Samples
for i=35:121
    disp(j)
    if j==13 || j==20 || j==66 || j==69
        AvgInd(i,:)=find(ismember(HILIC_NEG.All_Titles,['S' num2str(j) 'A']));
        if AvgInd(i,1)==(AvgInd(i-1,1))
            AvgInd(i,:)=find(ismember(HILIC_NEG.All_Titles,['S' num2str(j) 'B']));
            j=j+1;
        end
        
    elseif length(find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)])))==2
        AvgInd(i,1:2)=find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)]));
        j=j+1;
    elseif length(find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)])))==3
        AvgInd(i,:)=find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)]));
        j=j+1;
        
    else
        AvgInd(i,:)=0;
        j=j+1;
    end
end
r=find(sum(AvgInd,2)==0);
AvgInd(r,:)=[];

% Check it
 
 for i=1:length(AvgInd)
     if AvgInd(i,2)~=0 && AvgInd(i,3)~=0;
     A(i,1:3)=HILIC_NEG.All_Titles(AvgInd(i,:));
     elseif AvgInd(i,2)~=0;
     A(i,1:2)=HILIC_NEG.All_Titles(AvgInd(i,1:2));
     else
     A(i,1)=HILIC_NEG.All_Titles(AvgInd(i,1));   
     end
 end
 clear A r i j
 
% Step 2: Average along these rows 
for i=1:length(AvgInd)
    if AvgInd(i,2)~=0 && AvgInd(i,3)~=0
        HILIC_NEG.AvgNorm(i,:)=mean(HILIC_NEG.NormData(AvgInd(i,:),:),1);
    elseif AvgInd(i,2)~=0 && AvgInd(i,3)==0
        HILIC_NEG.AvgNorm(i,:)=mean(HILIC_NEG.NormData(AvgInd(i,1:2),:),1);
    else
        HILIC_NEG.AvgNorm(i,:)=HILIC_NEG.NormData(AvgInd(i,1),:);
    end
end
 
 HILIC_NEG.AvgInd=AvgInd;
 HILIC_NEG.AvgTitles=HILIC_NEG.All_Titles(HILIC_NEG.AvgInd(:,1));
HILIC_NEG.AvgGroups=HILIC_NEG.All_Groups(HILIC_NEG.AvgInd(:,1));

 clear A r i j AvgInd

  %% Average Replicates Positive
%Step 1: Create a matrix with the indeces corresponding to the replicates

%Step 1: Create a matrix with the indeces corresponding to the replicates

AvgInd=zeros(100,3);

for i=1:11 % Pooled
AvgInd(i,1:3)=[i 0 0];
end

HILIC_POS.All_Titles=strrep(HILIC_POS.All_Titles, '_1', '');
HILIC_POS.All_Titles=strrep(HILIC_POS.All_Titles, '_2', '');
HILIC_POS.All_Titles=strrep(HILIC_POS.All_Titles, '_3', '');
%%
j=1; %Samples
for i=12:100
    disp(j)
    if j==13 || j==20 || j==66 || j==69
         AvgInd(i,:)=find(ismember(HILIC_POS.All_Titles,['S' num2str(j) 'A']));    
        if AvgInd(i,1)==(AvgInd(i-1,1))
            AvgInd(i,:)=find(ismember(HILIC_POS.All_Titles,['S' num2str(j) 'B']));
            j=j+1;
        end
        
    elseif ~isempty(find(ismember(HILIC_POS.All_Titles,['S' num2str(j)])))
        AvgInd(i,:)=find(ismember(HILIC_POS.All_Titles,['S' num2str(j)]));
        j=j+1;
    else
        AvgInd(i,:)=0;
        j=j+1;
    end
end
r=find(sum(AvgInd,2)==0);
AvgInd(r,:)=[];

%CHECK

for i=1:length(AvgInd)
     if AvgInd(i,3)~=0 
     A(i,1:3)=HILIC_POS.All_Titles(AvgInd(i,:));
     else
     A(i,1)=HILIC_POS.All_Titles(AvgInd(i,1));   
     end
 end
 clear A r i j
 
% Step 2: Average along these rows 
 for i=1:length(AvgInd)
     if AvgInd(i,2)~=0
         HILIC_POS.AvgNorm(i,:)=mean(HILIC_POS.NormData(AvgInd(i,:),:),1); 
     else
         HILIC_POS.AvgNorm(i,:)=HILIC_POS.NormData(AvgInd(i,1),:); 
     end
 end
 
 HILIC_POS.AvgInd=AvgInd;
 clear A r i j AvgInd
 HILIC_POS.AvgTitles=HILIC_POS.All_Titles(HILIC_POS.AvgInd(:,1));
 HILIC_POS.AvgGroups=HILIC_POS.All_Groups(HILIC_POS.AvgInd(:,1));

 %% Boxplot to Visualize
normcheckSubplot(HILIC_NEG.NormData,HILIC_NEG.AvgNorm,HILIC_NEG.AvgTitles,'one','Normalization','two','Average Normalization');
normcheckSubplot(HILIC_POS.NormData,HILIC_POS.AvgNorm,HILIC_POS.AvgTitles,'one','Normalization','two','Average Normalization');

maboxplot(log((HILIC_POS.AvgNorm./repmat(median(HILIC_POS.AvgNorm),[size(HILIC_POS.AvgNorm,1),1]))'),'title','Averaged Progenesis Normalized Data','orientation','horizontal')
ylabel('log(Feature/Median Feature Value)')
xlabel('Sample')
set(gca,'XTickLabel',HILIC_POS.AvgTitles)

maboxplot(log((HILIC_NEG.AvgNorm./repmat(median(HILIC_NEG.AvgNorm),[size(HILIC_NEG.AvgNorm,1),1]))'),'title','Averaged Progenesis Normalized Data','orientation','horizontal')
ylabel('log(Feature/Median Feature Value)')
xlabel('Sample')
set(gca,'XTickLabel',HILIC_NEG.AvgTitles)

%% Remove Pooled, QC, Sample 77 
HILIC_NEG.AvgNorm_wQP=HILIC_NEG.AvgNorm(35:end,:);
HILIC_NEG.AvgTitles_wQP=HILIC_NEG.AvgTitles(35:end);
HILIC_NEG.AvgGroups_wQP=HILIC_NEG.AvgGroups(35:end);

%% Remove Pooled, QC, Sample 77 
HILIC_POS.AvgNorm_wQP=HILIC_POS.AvgNorm(12:end,:);
HILIC_POS.AvgTitles_wQP=HILIC_POS.AvgTitles(12:end);
HILIC_POS.AvgGroups_wQP=HILIC_POS.AvgGroups(12:end);

%% Create Y Vector Negative
HILIC_NEG.Y=zeros(length(HILIC_NEG.AvgGroups_wQP),1);
HILIC_NEG.Y(find(ismember(HILIC_NEG.AvgGroups_wQP,'No Recurr')))=2;  %2 is No Reccurance, 1 is Recurrance
HILIC_NEG.Y(find(ismember(HILIC_NEG.AvgGroups_wQP,'Recurr')))=1;

HILIC_NEG.ID=HILIC_NEG.DataFilt(3:end,[6:8 10 4]);
%% Get the significance of each Data set with an FDR correction
%HILIC NEG
norecurr=HILIC_NEG.AvgNorm_wQP(HILIC_NEG.Y==2,:); %no recur
recurr=HILIC_NEG.AvgNorm_wQP(HILIC_NEG.Y==1,:); %recur
norecurr_tit=HILIC_NEG.AvgTitles_wQP(HILIC_NEG.Y==2); %no recur
recurr_tit=HILIC_NEG.AvgTitles_wQP(HILIC_NEG.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

for i=1:length(norecurr)
    [~,HILIC_NEG.Sign.Pvalues(i)]=ttest2(norecurr(:,i),recurr(:,i),'Vartype','unequal');
end

[~, HILIC_NEG.AdjSign.Crit_pvalue, HILIC_NEG.AdjSign.Adjpvalues]=fdr_bh(HILIC_NEG.Sign.Pvalues); %if using FDR 
HILIC_NEG.Sign.Signindx=find(HILIC_NEG.Sign.Pvalues<0.05);
HILIC_NEG.Sign.SignInt=HILIC_NEG.AvgNorm_wQP(:,HILIC_NEG.Sign.Signindx);
HILIC_NEG.Sign.SignMRT=HILIC_NEG.MASS_RT(:,HILIC_NEG.Sign.Signindx)';
HILIC_NEG.Sign.SignID=HILIC_NEG.ID([1 HILIC_NEG.Sign.Signindx+1],:);

HILIC_NEG.AdjSign.AdjSignindx=find(HILIC_NEG.AdjSign.Adjpvalues<0.05);
HILIC_NEG.AdjSign.AdjSignInt=HILIC_NEG.AvgNorm_wQP(:,HILIC_NEG.AdjSign.AdjSignindx);
HILIC_NEG.AdjSign.AdjSignMRT=HILIC_NEG.MASS_RT(:,HILIC_NEG.AdjSign.AdjSignindx)';
HILIC_NEG.AdjSign.AdjSignID=HILIC_NEG.ID([1 HILIC_NEG.AdjSign.AdjSignindx+1],:);


HILIC_NEG.Sign.SignFoldChange=VolcanoPlot(norecurr,recurr,HILIC_NEG.Sign.Pvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',HILIC_NEG.MASS_RT(1,:),'Features2',HILIC_NEG.MASS_RT(2,:),'Groups',Groups);
HILIC_NEG.Sign.SignFoldChange.ID=HILIC_NEG.ID([1 HILIC_NEG.Sign.SignFoldChange.idx+1],:);

HILIC_NEG.AdjSign.AdjFoldChange=VolcanoPlot(norecurr,recurr,HILIC_NEG.AdjSign.Adjpvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',HILIC_NEG.MASS_RT(1,:),'Features2',HILIC_NEG.MASS_RT(2,:),'Groups',Groups);
HILIC_NEG.AdjSign.AdjFoldChange.ID=HILIC_NEG.ID([1 HILIC_NEG.AdjSign.AdjFoldChange.idx+1],:);

clear norecurr recurr recurr_tit norecurr_tit Groups i

%% Create Y Vector Positive
HILIC_POS.Y=zeros(length(HILIC_POS.AvgGroups_wQP),1);
HILIC_POS.Y(find(ismember(HILIC_POS.AvgGroups_wQP,'No Recurr')))=2;  %2 is No Reccurance, 1 is Recurrance
HILIC_POS.Y(find(ismember(HILIC_POS.AvgGroups_wQP,'Recurr')))=1;

HILIC_POS.ID=HILIC_POS.DataFilt(3:end,[6:8 10 4]);
%% Get the significance of each Data set with an FDR correction
%HILIC POS
norecurr=HILIC_POS.AvgNorm_wQP(HILIC_POS.Y==2,:); %no recur
recurr=HILIC_POS.AvgNorm_wQP(HILIC_POS.Y==1,:); %recur
norecurr_tit=HILIC_POS.AvgTitles_wQP(HILIC_POS.Y==2); %no recur
recurr_tit=HILIC_POS.AvgTitles_wQP(HILIC_POS.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

for i=1:length(norecurr)
    [~,HILIC_POS.Sign.Pvalues(i)]=ttest2(norecurr(:,i),recurr(:,i),'Vartype','unequal');
end

[~, HILIC_POS.AdjSign.Crit_pvalue, HILIC_POS.AdjSign.Adjpvalues]=fdr_bh(HILIC_POS.Sign.Pvalues); %if using FDR 
HILIC_POS.Sign.Signindx=find(HILIC_POS.Sign.Pvalues<0.05);
HILIC_POS.Sign.SignInt=HILIC_POS.AvgNorm_wQP(:,HILIC_POS.Sign.Signindx);
HILIC_POS.Sign.SignMRT=HILIC_POS.MASS_RT(:,HILIC_POS.Sign.Signindx)';
HILIC_POS.Sign.SignID=HILIC_POS.ID([1 HILIC_POS.Sign.Signindx+1],:);

HILIC_POS.AdjSign.AdjSignindx=find(HILIC_POS.AdjSign.Adjpvalues<0.05);
HILIC_POS.AdjSign.AdjSignInt=HILIC_POS.AvgNorm_wQP(:,HILIC_POS.AdjSign.AdjSignindx);
HILIC_POS.AdjSign.AdjSignMRT=HILIC_POS.MASS_RT(:,HILIC_POS.AdjSign.AdjSignindx)';
HILIC_POS.AdjSign.AdjSignID=HILIC_POS.ID([1 HILIC_POS.AdjSign.AdjSignindx+1],:);


HILIC_POS.Sign.SignFoldChange=VolcanoPlot(norecurr,recurr,HILIC_POS.Sign.Pvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',HILIC_POS.MASS_RT(1,:),'Features2',HILIC_POS.MASS_RT(2,:),'Groups',Groups);
HILIC_POS.Sign.SignFoldChange.ID=HILIC_POS.ID([1 HILIC_POS.Sign.SignFoldChange.idx+1],:);

HILIC_POS.AdjSign.AdjFoldChange=VolcanoPlot(norecurr,recurr,HILIC_POS.AdjSign.Adjpvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',HILIC_POS.MASS_RT(1,:),'Features2',HILIC_POS.MASS_RT(2,:),'Groups',Groups);
HILIC_POS.AdjSign.AdjFoldChange.ID=HILIC_POS.ID([1 HILIC_POS.AdjSign.AdjFoldChange.idx+1],:);

clear norecurr reccur norecurr_tit reccur_tit Groups i

%% Fold Change
norecurr=HILIC_POS.AvgNorm(HILIC_POS.Y==2,:); %no recur
recurr=HILIC_POS.AvgNorm(HILIC_POS.Y==1,:); %recur
norecurr_tit=HILIC_POS.AvgTitles_wQP(HILIC_POS.Y==2); %no recur
recurr_tit=HILIC_POS.AvgTitles_wQP(HILIC_POS.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

HILIC_POS.FoldChange=mean(norecurr)./mean(recurr);
HILIC_POS.Log2FC=log2(HILIC_POS.FoldChange);