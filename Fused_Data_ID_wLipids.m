%% Create a Y1 and Y2 vector that describes the classes
% NMR.All_Titles=strrep(NMR.All_Titles,'''','');
All_Titles=NMR.All_Titles;
Fused_Data.Y=zeros(68,1);
Fused_Data.Y(find(ismember(All_Titles,'S1')))=2;  %2 is No Recurrence, 1 is Recurrence
Fused_Data.Y(find(ismember(All_Titles,'S2')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S3')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S4')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S5')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S6')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S7')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S8')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S9')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S12')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S13')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S14')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S16')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S17')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S18')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S19')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S21')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S22')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S23')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S24')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S26')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S27')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S28')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S29')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S30')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S31')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S32')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S33')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S36')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S38')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S41')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S42')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S44')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S46')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S47')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S66')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S67')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S70')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S72')))=2;
Fused_Data.Y(find(ismember(All_Titles,'S77')))=2;
Fused_Data.Y(Fused_Data.Y==0)=1;

%% Combine the 68 NMR samples and features with 68 of the RP_POS and NEG
% and HILIC
% samples

for i=1:length(NMR.All_Titles)  %Use Sample A for the double
    if ismember(NMR.All_Titles(i),'S13')
        [~,idx]=ismember(HILIC_POS.AvgTitles_wQP,'S13A');
        indx(2,i)=find(idx==1);
    elseif ismember(NMR.All_Titles(i),'S20')
        [~,idx]=ismember(HILIC_POS.AvgTitles_wQP,'S20A');
        indx(2,i)=find(idx==1);
    elseif ismember(NMR.All_Titles(i),'S66')
        [~,idx]=ismember(HILIC_POS.AvgTitles_wQP,'S66A');
        indx(2,i)=find(idx==1);
    elseif ismember(NMR.All_Titles(i),'S69')
        [~,idx]=ismember(HILIC_POS.AvgTitles_wQP,'S69A');
        indx(2,i)=find(idx==1);
    else
        [~,idx]=ismember(NMR.All_Titles(i),HILIC_POS.AvgTitles_wQP);
        indx(2,i)=idx;
    end
    
end

clear i idx

%% Combine the 68 NMR samples and features with 68 of the RP_POS and NEG
% and RP
% samples

for i=1:length(NMR.All_Titles)  %Use Sample A for the double
        [~,idx]=ismember(NMR.All_Titles(i),RP_POS.AvgTitles_wQP);
        indx(2,i)=idx;   
end

clear i idx

%% The samples that has not been in the NMR set
RP=1:79;
[idx]=setdiff(RP,indx(1,:));

clear RP HILIC

Fused_Data.UnFused_Data.RP=RP_POS.AvgTitles_wQP(idx);

Fused_Data.UnFused_Data.RP_POS=RP_POS.AvgNorm_wQP(idx,:);
Fused_Data.UnFused_Data.RP_NEG=RP_NEG.AvgNorm_wQP(idx,:);

Fused_Data.UnFused_Data.RP_POS_MRT=RP_POS.MASS_RT;
Fused_Data.UnFused_Data.RP_NEG_MRT=RP_NEG.MASS_RT;

Fused_Data.UnFused_Data.RP_POS_Class=RP_POS.AvgGroups_wQP(idx);
Fused_Data.UnFused_Data.RP_NEG_Class=RP_NEG.AvgGroups_wQP(idx);

Fused_Data.UnFused_Data.HILIC=HILIC_POS.AvgTitles_wQP(idx2);

Fused_Data.UnFused_Data.HILIC_POS=HILIC_POS.AvgNorm_wQP(idx2,:);
Fused_Data.UnFused_Data.HILIC_NEG=HILIC_NEG.AvgNorm_wQP(idx2,:);

Fused_Data.UnFused_Data.HILIC_POS_MRT=HILIC_POS.MASS_RT;
Fused_Data.UnFused_Data.HILIC_NEG_MRT=HILIC_NEG.MASS_RT;

Fused_Data.UnFused_Data.HILIC_POS_Class=HILIC_POS.AvgGroups_wQP(idx2);
Fused_Data.UnFused_Data.HILIC_NEG_Class=HILIC_NEG.AvgGroups_wQP(idx2);

%% Make a new structure called Fused_Data
Fused_Data.Fused_Titles_HILIC=HILIC_POS.AvgTitles_wQP(indx(2,:));
Fused_Data.Fused_Titles_RP=RP_POS.AvgTitles_wQP(indx(2,:));
Fused_Data.Fused_Titles_NMR=NMR.All_Titles';

Fused_Data.HILIC_POS.Intensity=HILIC_POS.AvgNorm_wQP(indx(2,:),:);
Fused_Data.HILIC_NEG.Intensity=HILIC_NEG.AvgNorm_wQP(indx(2,:),:);
Fused_Data.RP_POS.Intensity=RP_POS.AvgNorm_wQP(indx(2,:),:);
Fused_Data.RP_NEG.Intensity=RP_NEG.AvgNorm_wQP(indx(2,:),:);
Fused_Data.NMR.Intensity=NMR.NormData;
Fused_Data.NMR_sum.Intensity=NMR.NormData_sum;

Fused_Data.HILIC_POS.MRT=HILIC_POS.MASS_RT;
Fused_Data.HILIC_NEG.MRT=HILIC_NEG.MASS_RT;
Fused_Data.RP_POS.MRT=RP_POS.MASS_RT;
Fused_Data.RP_NEG.MRT=RP_NEG.MASS_RT;
Fused_Data.NMR.ppm=NMR.ppm;
Fused_Data.NMR_sum.ppm=NMR.ppm_sum;

Fused_Data.HILIC_POS.ID=HILIC_POS.ID(2:end,2);
Fused_Data.HILIC_NEG.ID=HILIC_NEG.ID(2:end,1);
Fused_Data.RP_POS.ID=RP_POS.ID(2:end,3);
Fused_Data.RP_NEG.ID=RP_NEG.ID(2:end,3);
Fused_Data.NMR.ID=NMR.ID;
Fused_Data.NMR_sum.ID=NMR.ID_sum;

Fused_Data.HILIC_POS.Pvalues=HILIC_POS.Sign.Pvalues;
Fused_Data.HILIC_NEG.Pvalues=HILIC_NEG.Sign.Pvalues;
Fused_Data.RP_POS.Pvalues=RP_POS.Sign.Pvalues;
Fused_Data.RP_NEG.Pvalues=RP_NEG.Sign.Pvalues;
Fused_Data.NMR.Pvalues=NMR.Sign.Pvalues;
Fused_Data.NMR_sum.Pvalues=NMR.Sign_sum.Pvalues;

Fused_Data.HILIC_POS.AdjPvalues=HILIC_POS.AdjSign.Adjpvalues;
Fused_Data.HILIC_NEG.AdjPvalues=HILIC_NEG.AdjSign.Adjpvalues;
Fused_Data.RP_POS.AdjPvalues=RP_POS.AdjSign.Adjpvalues;
Fused_Data.RP_NEG.AdjPvalues=RP_NEG.AdjSign.Adjpvalues;
Fused_Data.NMR.AdjPvalues=NMR.AdjSign.Adjpvalues;
Fused_Data.NMR_sum.AdjPvalues=NMR.AdjSign_sum.Adjpvalues;


clear idx idx2 RP HILIC indx

%% Create a Fused Matrix
Fused_Data.Fused_Data.All_Intensity=[Fused_Data.HILIC_POS.Intensity...
    Fused_Data.HILIC_NEG.Intensity Fused_Data.RP_POS.Intensity...
    Fused_Data.RP_NEG.Intensity Fused_Data.NMR.Intensity];
Fused_Data.Fused_Data.All_Points=[Fused_Data.HILIC_POS.MRT(1,:) Fused_Data.HILIC_NEG.MRT(1,:)... 
    Fused_Data.RP_POS.MRT(1,:) Fused_Data.RP_NEG.MRT(1,:)...
    Fused_Data.NMR.ppm; Fused_Data.HILIC_POS.MRT(2,:)...
    Fused_Data.HILIC_NEG.MRT(2,:) Fused_Data.RP_POS.MRT(2,:)...
    Fused_Data.RP_NEG.MRT(2,:) Fused_Data.NMR.ppm];
Fused_Data.Fused_Data.Fused_Classes=NMR.All_Groups;
Fused_Data.Fused_Data.Fused_Titles=NMR.All_Titles';

Fused_Data.Fused_Data.Source={};
Fused_Data.Fused_Data.Source(end+1:end+length(Fused_Data.HILIC_POS.MRT))={'HILIC_POS'};
Fused_Data.Fused_Data.Source(end+1:end+length(Fused_Data.HILIC_NEG.MRT))={'HILIC_NEG'};
Fused_Data.Fused_Data.Source(end+1:end+length(Fused_Data.RP_POS.MRT))={'RP_POS'};
Fused_Data.Fused_Data.Source(end+1:end+length(Fused_Data.RP_NEG.MRT))={'RP_NEG'};
Fused_Data.Fused_Data.Source(end+1:end+length(Fused_Data.NMR.ppm))={'NMR'};

Fused_Data.Fused_Data.ID={};
Fused_Data.Fused_Data.ID(end+1:end+length(Fused_Data.HILIC_POS.ID))=Fused_Data.HILIC_POS.ID;
Fused_Data.Fused_Data.ID(end+1:end+length(Fused_Data.HILIC_NEG.ID))=Fused_Data.HILIC_NEG.ID;
Fused_Data.Fused_Data.ID(end+1:end+length(Fused_Data.RP_POS.ID))=Fused_Data.RP_POS.ID;
Fused_Data.Fused_Data.ID(end+1:end+length(Fused_Data.RP_NEG.ID))=Fused_Data.RP_NEG.ID;
Fused_Data.Fused_Data.ID(end+1:end+length(Fused_Data.NMR.ID))=Fused_Data.NMR.ID;

Fused_Data.Fused_Data.Pvalues=[];
Fused_Data.Fused_Data.Pvalues(end+1:end+length(Fused_Data.HILIC_POS.Pvalues))=Fused_Data.HILIC_POS.Pvalues;
Fused_Data.Fused_Data.Pvalues(end+1:end+length(Fused_Data.HILIC_NEG.Pvalues))=Fused_Data.HILIC_NEG.Pvalues;
Fused_Data.Fused_Data.Pvalues(end+1:end+length(Fused_Data.RP_POS.Pvalues))=Fused_Data.RP_POS.Pvalues;
Fused_Data.Fused_Data.Pvalues(end+1:end+length(Fused_Data.RP_NEG.Pvalues))=Fused_Data.RP_NEG.Pvalues;
Fused_Data.Fused_Data.Pvalues(end+1:end+length(Fused_Data.NMR.Pvalues))=Fused_Data.NMR.Pvalues;

Fused_Data.Fused_Data.AdjPvalues=[];
Fused_Data.Fused_Data.AdjPvalues(end+1:end+length(Fused_Data.HILIC_POS.AdjPvalues))=Fused_Data.HILIC_POS.AdjPvalues;
Fused_Data.Fused_Data.AdjPvalues(end+1:end+length(Fused_Data.HILIC_NEG.AdjPvalues))=Fused_Data.HILIC_NEG.AdjPvalues;
Fused_Data.Fused_Data.AdjPvalues(end+1:end+length(Fused_Data.RP_POS.AdjPvalues))=Fused_Data.RP_POS.AdjPvalues;
Fused_Data.Fused_Data.AdjPvalues(end+1:end+length(Fused_Data.RP_NEG.AdjPvalues))=Fused_Data.RP_NEG.AdjPvalues;
Fused_Data.Fused_Data.AdjPvalues(end+1:end+length(Fused_Data.NMR.AdjPvalues))=Fused_Data.NMR.AdjPvalues;

%% Get the Fold Change

norecurr=Fused_Data.Fused_Data.All_Intensity(Fused_Data.Y==2,:); %no recur
recurr=Fused_Data.Fused_Data.All_Intensity(Fused_Data.Y==1,:); %recur
norecurr_tit=Fused_Data.Fused_Data.Fused_Titles(Fused_Data.Y==2); %no recur
recurr_tit=Fused_Data.Fused_Data.Fused_Titles(Fused_Data.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

Fused_Data.Fused_Data.FoldChange=mean(norecurr)./mean(recurr);
Fused_Data.Fused_Data.Log2FC=log2(Fused_Data.Fused_Data.FoldChange);

%% Determine correlations

[Fused_Data.Fused_Data.RHO,Fused_Data.Fused_Data.PVAL]=corr(Fused_Data.Fused_Data.All_Intensity);
imagesc(Fused_Data.Fused_Data.RHO);
colormap(redblue)

%% Create a Fused Matrix with small nmr set
Fused_Data.Fused_Data_sum.All_Intensity=[Fused_Data.HILIC_POS.Intensity...
    Fused_Data.HILIC_NEG.Intensity Fused_Data.RP_POS.Intensity...
    Fused_Data.RP_NEG.Intensity Fused_Data.NMR_sum.Intensity];
Fused_Data.Fused_Data_sum.All_Points=[Fused_Data.HILIC_POS.MRT(1,:) Fused_Data.HILIC_NEG.MRT(1,:)... 
    Fused_Data.RP_POS.MRT(1,:) Fused_Data.RP_NEG.MRT(1,:)...
    Fused_Data.NMR_sum.ppm; Fused_Data.HILIC_POS.MRT(2,:)...
    Fused_Data.HILIC_NEG.MRT(2,:) Fused_Data.RP_POS.MRT(2,:)...
    Fused_Data.RP_NEG.MRT(2,:) Fused_Data.NMR_sum.ppm];
Fused_Data.Fused_Data_sum.Fused_Classes=NMR.All_Groups;
Fused_Data.Fused_Data_sum.Fused_Titles=NMR.All_Titles';

Fused_Data.Fused_Data_sum.Source={};
Fused_Data.Fused_Data_sum.Source(end+1:end+length(Fused_Data.HILIC_POS.MRT))={'HILIC_POS'};
Fused_Data.Fused_Data_sum.Source(end+1:end+length(Fused_Data.HILIC_NEG.MRT))={'HILIC_NEG'};
Fused_Data.Fused_Data_sum.Source(end+1:end+length(Fused_Data.RP_POS.MRT))={'RP_POS'};
Fused_Data.Fused_Data_sum.Source(end+1:end+length(Fused_Data.RP_NEG.MRT))={'RP_NEG'};
Fused_Data.Fused_Data_sum.Source(end+1:end+length(Fused_Data.NMR_sum.ppm))={'NMR'};

Fused_Data.Fused_Data_sum.ID={};
Fused_Data.Fused_Data_sum.ID(end+1:end+length(Fused_Data.HILIC_POS.ID))=Fused_Data.HILIC_POS.ID;
Fused_Data.Fused_Data_sum.ID(end+1:end+length(Fused_Data.HILIC_NEG.ID))=Fused_Data.HILIC_NEG.ID;
Fused_Data.Fused_Data_sum.ID(end+1:end+length(Fused_Data.RP_POS.ID))=Fused_Data.RP_POS.ID;
Fused_Data.Fused_Data_sum.ID(end+1:end+length(Fused_Data.RP_NEG.ID))=Fused_Data.RP_NEG.ID;
Fused_Data.Fused_Data_sum.ID(end+1:end+length(Fused_Data.NMR_sum.ID))=Fused_Data.NMR_sum.ID;

Fused_Data.Fused_Data_sum.Pvalues=[];
Fused_Data.Fused_Data_sum.Pvalues(end+1:end+length(Fused_Data.HILIC_POS.Pvalues))=Fused_Data.HILIC_POS.Pvalues;
Fused_Data.Fused_Data_sum.Pvalues(end+1:end+length(Fused_Data.HILIC_NEG.Pvalues))=Fused_Data.HILIC_NEG.Pvalues;
Fused_Data.Fused_Data_sum.Pvalues(end+1:end+length(Fused_Data.RP_POS.Pvalues))=Fused_Data.RP_POS.Pvalues;
Fused_Data.Fused_Data_sum.Pvalues(end+1:end+length(Fused_Data.RP_NEG.Pvalues))=Fused_Data.RP_NEG.Pvalues;
Fused_Data.Fused_Data_sum.Pvalues(end+1:end+length(Fused_Data.NMR_sum.Pvalues))=Fused_Data.NMR_sum.Pvalues;

Fused_Data.Fused_Data_sum.AdjPvalues=[];
Fused_Data.Fused_Data_sum.AdjPvalues(end+1:end+length(Fused_Data.HILIC_POS.AdjPvalues))=Fused_Data.HILIC_POS.AdjPvalues;
Fused_Data.Fused_Data_sum.AdjPvalues(end+1:end+length(Fused_Data.HILIC_NEG.AdjPvalues))=Fused_Data.HILIC_NEG.AdjPvalues;
Fused_Data.Fused_Data_sum.AdjPvalues(end+1:end+length(Fused_Data.RP_POS.AdjPvalues))=Fused_Data.RP_POS.AdjPvalues;
Fused_Data.Fused_Data_sum.AdjPvalues(end+1:end+length(Fused_Data.RP_NEG.AdjPvalues))=Fused_Data.RP_NEG.AdjPvalues;
Fused_Data.Fused_Data_sum.AdjPvalues(end+1:end+length(Fused_Data.NMR_sum.AdjPvalues))=Fused_Data.NMR_sum.AdjPvalues;

%% Get the Fold Change

norecurr=Fused_Data.Fused_Data_sum.All_Intensity(Fused_Data.Y==2,:); %no recur
recurr=Fused_Data.Fused_Data_sum.All_Intensity(Fused_Data.Y==1,:); %recur
norecurr_tit=Fused_Data.Fused_Data_sum.Fused_Titles(Fused_Data.Y==2); %no recur
recurr_tit=Fused_Data.Fused_Data_sum.Fused_Titles(Fused_Data.Y==1); %recur
Groups{1}='No Recurrence';
Groups{2}='Recurrence';

Fused_Data.Fused_Data_sum.FoldChange=mean(norecurr)./mean(recurr);
Fused_Data.Fused_Data_sum.Log2FC=log2(Fused_Data.Fused_Data_sum.FoldChange);
%% Determine correlations

[Fused_Data.Fused_Data_sum.RHO,Fused_Data.Fused_Data_sum.PVAL]=corr(Fused_Data.Fused_Data_sum.All_Intensity,'type','Spearman');
imagesc(Fused_Data.Fused_Data_sum.RHO);

colormap(redblue)
%% Determine correlations
%Remove column 81 cause isnan
Fused_Data.Fused_Data_sum.All_Intensity_81=Fused_Data.Fused_Data_sum.All_Intensity(:,[1:80 82:end]);
Fused_Data.Fused_Data_sum.ID_81=Fused_Data.Fused_Data_sum.ID([1:80 82:end]);
Fused_Data.Fused_Data_sum.Source_81=Fused_Data.Fused_Data_sum.Source([1:80 82:end]);

[Fused_Data.Fused_Data_sum.RHO,Fused_Data.Fused_Data_sum.PVAL]=corr(Fused_Data.Fused_Data_sum.All_Intensity_81);
imagesc(Fused_Data.Fused_Data_sum.RHO);

colormap(redblue)
