%% Normalized Data imported as NMRID
% Create Structure
NMR.Data=NMRID;

%% Create a grouping list (Blank, Pooled, QC preremoved)

for i=3:size(NMR.Data,2)
    if ismember(NMR.Data(1,i),'''No Reccurance''')
        NMR.All_Groups{i}='No Recurr';
    elseif ismember(NMR.Data(1,i),'''Reccurance''')
        NMR.All_Groups{i}='Recurr';
    end
end

NMR.All_Groups(1:2)=[];

NMR.ppm=cell2mat(NMR.Data(3:end,1))';
NMR.ID=NMR.Data(3:end,2);
NMR.All_Titles=NMR.Data(2,3:end);

NMR.NormData=cell2mat(NMR.Data(3:end,3:end))';


%% Create Y Vector Positive
NMR.Y=zeros(length(NMR.All_Groups),1);
NMR.Y(find(ismember(NMR.All_Groups,'No Recurr')))=2;  %2 is No Reccurance, 1 is Recurrance
NMR.Y(find(ismember(NMR.All_Groups,'Recurr')))=1;

%% Get the significance of each Data set with an FDR correction
%HILIC POS
norecurr=NMR.NormData(NMR.Y==2,:); %no recur
recurr=NMR.NormData(NMR.Y==1,:); %recur
norecurr_tit=NMR.All_Titles(NMR.Y==2); %no recur
recurr_tit=NMR.All_Titles(NMR.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

for i=1:length(norecurr)
    [~,NMR.Sign.Pvalues(i)]=ttest2(norecurr(:,i),recurr(:,i),'Vartype','unequal');
end

[~, NMR.AdjSign.Crit_pvalue, NMR.AdjSign.Adjpvalues]=fdr_bh(NMR.Sign.Pvalues); %if using FDR 
NMR.Sign.Signindx=find(NMR.Sign.Pvalues<0.05);
NMR.Sign.SignInt=NMR.NormData(:,NMR.Sign.Signindx);
NMR.Sign.SignMRT=NMR.ppm(:,NMR.Sign.Signindx)';
NMR.Sign.SignID=NMR.ID(NMR.Sign.Signindx,:);

NMR.AdjSign.AdjSignindx=find(NMR.AdjSign.Adjpvalues<0.05);
NMR.AdjSign.AdjSignInt=NMR.NormData(:,NMR.AdjSign.AdjSignindx);
NMR.AdjSign.AdjSignMRT=NMR.ppm(:,NMR.AdjSign.AdjSignindx)';
NMR.AdjSign.AdjSignID=NMR.ID(NMR.AdjSign.AdjSignindx,:);


NMR.Sign.SignFoldChange=VolcanoPlot(norecurr,recurr,NMR.Sign.Pvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',NMR.ppm(1,:),'Groups',Groups);
NMR.Sign.SignFoldChange.ID=NMR.ID(NMR.Sign.SignFoldChange.idx,:);

NMR.AdjSign.AdjFoldChange=VolcanoPlot(norecurr,recurr,NMR.AdjSign.Adjpvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',NMR.ppm(1,:),'Groups',Groups);
NMR.AdjSign.AdjFoldChange.ID=NMR.ID(NMR.AdjSign.AdjFoldChange.idx,:);

clear norecurr recurr norecurr_tit recurr_tit Groups i

%% Group like Metabolites
[~,~,p]=unique(NMR.ID);
Y=logical(diff(p));
A=find(Y==1);

for i=1:length(A)
    if i==1
    NMR.NormData_sum(:,i)=sum(NMR.NormData(:,1:A(i)),2);
    NMR.ID_sum(i)=NMR.ID(A(i));
    NMR.ppm_sum(i)=median(NMR.ppm(1:A(i)));
    else
    NMR.NormData_sum(:,i)=sum(NMR.NormData(:,A(i-1):A(i)),2); 
    NMR.ID_sum(i)=NMR.ID(A(i));
    NMR.ppm_sum(i)=median(NMR.ppm(A(i-1):A(i)));
    end
end
%% Get the significance of each Data set with an FDR correction
%HILIC POS
norecurr=NMR.NormData_sum(NMR.Y==2,:); %no recur
recurr=NMR.NormData_sum(NMR.Y==1,:); %recur
norecurr_tit=NMR.All_Titles(NMR.Y==2); %no recur
recurr_tit=NMR.All_Titles(NMR.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

for i=1:length(norecurr)
    [~,NMR.Sign_sum.Pvalues(i)]=ttest2(norecurr(:,i),recurr(:,i),'Vartype','unequal');
end

[~, NMR.AdjSign_sum.Crit_pvalue, NMR.AdjSign_sum.Adjpvalues]=fdr_bh(NMR.Sign_sum.Pvalues); %if using FDR 
NMR.Sign_sum.Signindx=find(NMR.Sign_sum.Pvalues<0.05);
NMR.Sign_sum.SignInt=NMR.NormData_sum(:,NMR.Sign_sum.Signindx);
NMR.Sign_sum.SignMRT=NMR.ppm_sum(:,NMR.Sign_sum.Signindx)';
NMR.Sign_sum.SignID=NMR.ID_sum(NMR.Sign_sum.Signindx);

NMR.AdjSign_sum.AdjSignindx=find(NMR.AdjSign_sum.Adjpvalues<0.05);
NMR.AdjSign_sum.AdjSignInt=NMR.NormData_sum(:,NMR.AdjSign_sum.AdjSignindx);
NMR.AdjSign_sum.AdjSignMRT=NMR.ppm_sum(:,NMR.AdjSign_sum.AdjSignindx)';
NMR.AdjSign_sum.AdjSignID=NMR.ID_sum(NMR.AdjSign_sum.AdjSignindx);


NMR.Sign_sum.SignFoldChange=VolcanoPlot(norecurr,recurr,NMR.Sign_sum.Pvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',NMR.ppm_sum(1,:),'Groups',Groups);
NMR.Sign_sum.SignFoldChange.ID=NMR.ID_sum(NMR.Sign_sum.SignFoldChange.idx);

NMR.AdjSign_sum.AdjFoldChange=VolcanoPlot(norecurr,recurr,NMR.AdjSign_sum.Adjpvalues,'NamesX',norecurr_tit','NamesY',recurr_tit','Features1',NMR.ppm_sum(1,:),'Groups',Groups);
NMR.AdjSign_sum.AdjFoldChange.ID=NMR.ID_sum(NMR.AdjSign_sum.AdjFoldChange.idx);

clear norecurr recurr norecurr_tit recurr_tit Groups i A B ans p N Y
