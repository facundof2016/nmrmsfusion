%% Positive mode 
RP_POS.Data=RPPositive;
%% Fix Isotopes and low peak width
j=1;
for i=1:length(RP_POS.Data)
    if isnumeric(RP_POS.Data{i,4})
        idx(j)=i;
        j=j+1;
    elseif i>3 && ismember(RP_POS.Data{i,5}(1),'1')
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans 
RP_POS.IsotopeIssues=[RP_POS.Data(1:3,:);RP_POS.Data(idx,:)];
RP_POS.Data(idx,:)=[];
clear idx

j=1;
for i=1:length(RP_POS.Data)
    if i>3 && RP_POS.Data{i,3}<0.08
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans
RP_POS.LowPeakWidth=[RP_POS.Data(1:3,:);RP_POS.Data(idx,:)];
RP_POS.Data(idx,:)=[];
clear idx

%%
RP_POS.All_Groups=RP_POS.Data(2,5:185);
for i=2:length(RP_POS.All_Groups)
    if i<6
        RP_POS.All_Groups{i}='Blank';
    elseif i>=6 && i<12
        RP_POS.All_Groups{i}='Pooled';
    elseif i>=12 && i<23
        RP_POS.All_Groups{i}='QC';
    elseif i>=23 && i<102
        RP_POS.All_Groups{i}='No Recurrence';
    elseif i>=102
        RP_POS.All_Groups{i}='Recurrence';
    end
end

RP_POS.MASS_RT=cell2mat(RP_POS.Data(4:end,[1 2]))';
RP_POS.All_Titles=RP_POS.Data(3,5:185);
RP_POS.All_Titles(find(ismember(RP_POS.All_Groups,'Blank')))=[];

RP_POS.NormData=cell2mat(RP_POS.Data(4:end,5:185))';
RP_POS.NormData(find(ismember(RP_POS.All_Groups,'Blank')),:)=[];

RP_POS.All_Groups(find(ismember(RP_POS.All_Groups,'Blank')))=[];

%% Negative mode 
RP_NEG.Data=RPNegative;

%% Fix Isotopes and low peak width
j=1;
for i=1:length(RP_NEG.Data)
    if isnumeric(RP_NEG.Data{i,4})
        idx(j)=i;
        j=j+1;
    elseif i>3 && ismember(RP_NEG.Data{i,5}(1),'1')
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans 
RP_NEG.IsotopeIssues=[RP_NEG.Data(1:3,:);RP_NEG.Data(idx,:)];
RP_NEG.Data(idx,:)=[];
clear idx

j=1;
for i=1:length(RP_NEG.Data)
    if i>3 &&RP_NEG.Data{i,3}<0.08
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans
RP_NEG.LowPeakWidth=[RP_NEG.Data(1:3,:);RP_NEG.Data(idx,:)];
RP_NEG.Data(idx,:)=[];
clear idx
%%
RP_NEG.All_Groups=RP_NEG.Data(2,5:188);
for i=1:length(RP_NEG.All_Groups)
    if i<12
        RP_NEG.All_Groups{i}='Blank';
    elseif i>=12 && i<16
        RP_NEG.All_Groups{i}='Pooled';
    elseif i>=16 && i<25
        RP_NEG.All_Groups{i}='QC';
    elseif i>=25 && i<105
        RP_NEG.All_Groups{i}='No Recurr';
    elseif i>=105
        RP_NEG.All_Groups{i}='Recurr';
    end
end

RP_NEG.MASS_RT=cell2mat(RP_NEG.Data(4:end,[1 3]))';
RP_NEG.All_Titles=RP_NEG.Data(3,5:188);
RP_NEG.All_Titles(find(ismember(RP_NEG.All_Groups,'Blank')))=[];

RP_NEG.NormData=cell2mat(RP_NEG.Data(4:end,5:188))';
RP_NEG.NormData(find(ismember(RP_NEG.All_Groups,'Blank')),:)=[];

RP_NEG.All_Groups(find(ismember(RP_NEG.All_Groups,'Blank')))=[];

%% Remove first few seconds from each of the MS methods
idx1=find(RP_POS.MASS_RT(2,:)>0.65);

idx2=find(RP_NEG.MASS_RT(2,:)>0.65);

RP_POS.MASS_RT_cut=RP_POS.MASS_RT(:,idx1);
RP_NEG.MASS_RT_cut=RP_NEG.MASS_RT(:,idx2);

RP_POS.NormData_cut=RP_POS.NormData(:,idx1);
RP_NEG.NormData_cut=RP_NEG.NormData(:,idx2);
clear idx1 idx2

idx1=find(RP_POS.MASS_RT(2,:)<=0.65);
idx2=find(RP_NEG.MASS_RT(2,:)<=0.65);
RP_POS.Cut.MASS_RT=RP_POS.MASS_RT(:,idx1);
RP_NEG.Cut.MASS_RT=RP_NEG.MASS_RT(:,idx2);

RP_POS.Cut.NormData=RP_POS.NormData(:,idx1);
RP_NEG.Cut.NormData=RP_NEG.NormData(:,idx2);

clear idx2 idx1

%% Average Replicates Negative
%Step 1: Create a matrix with the indeces corresponding to the replicates

AvgInd_cut=zeros(length(RP_NEG.All_Titles),2);

for i=1:13 % Pooled and QC
AvgInd_cut(i,1:2)=[i 0];
end

RP_NEG.All_Titles=strrep(RP_NEG.All_Titles, '_1', '');
RP_NEG.All_Titles=strrep(RP_NEG.All_Titles, '_2', '');
%%
j=1; %Samples
for i=14:173
    disp(j)
    if ~isempty(find(ismember(RP_NEG.All_Titles,['S' num2str(j)])));
        AvgInd_cut(i,:)=find(ismember(RP_NEG.All_Titles,['S' num2str(j)]));
        j=j+1;
   
    else
        AvgInd_cut(i,:)=0;
        j=j+1;
    end
end
[r,~]=find(sum(AvgInd_cut,2)==0);
AvgInd_cut(r,:)=[];

% Step 2: Average along these rows 
  for i=1:length(AvgInd_cut)
     if AvgInd_cut(i,2)~=0
         RP_NEG.AvgNorm_cut(i,:)=mean(RP_NEG.NormData_cut(AvgInd_cut(i,:),:),1); 
         RP_NEG.AvgTitles(i)=RP_NEG.All_Titles(AvgInd_cut(i,1));
     else
         RP_NEG.AvgNorm_cut(i,:)=RP_NEG.NormData_cut(AvgInd_cut(i,1),:); 
         RP_NEG.AvgTitles(i)=RP_NEG.All_Titles(AvgInd_cut(i,1));
     end
 end
 
 RP_NEG.AvgInd_cut=AvgInd_cut;
 RP_NEG.AvgGroups=RP_NEG.All_Groups(RP_NEG.AvgInd_cut(:,1));
 clear A r i j AvgInd_cut

   %% Average Replicates Positive
%Step 1: Create a matrix with the indeces corresponding to the replicates

AvgInd_cut=zeros(length(RP_POS.All_Titles),2);

for i=1:16 % Pooled and QC
    AvgInd_cut(i,1:2)=[i 0];
end

RP_POS.All_Titles=strrep(RP_POS.All_Titles, '_1', '');
RP_POS.All_Titles=strrep(RP_POS.All_Titles, '_2', '');

%%
j=1; %Samples
for i=17:176
    %disp(j)
    if ~isempty(find(ismember(RP_POS.All_Titles,['S' num2str(j)])));
        AvgInd_cut(i,:)=find(ismember(RP_POS.All_Titles,['S' num2str(j)]));
    else
         AvgInd_cut(i,1:2)=0;         
    end
    j=j+1;
end
[r,~]=find(sum(AvgInd_cut,2)==0);
AvgInd_cut(r,:)=[];


% Step 2: Average along these rows 
 for i=1:length(AvgInd_cut)
     if AvgInd_cut(i,2)~=0
         RP_POS.AvgNorm_cut(i,:)=mean(RP_POS.NormData_cut(AvgInd_cut(i,:),:),1); 
         RP_POS.AvgTitles(i)=RP_POS.All_Titles(AvgInd_cut(i,1));
     else
         RP_POS.AvgNorm_cut(i,:)=RP_POS.NormData_cut(AvgInd_cut(i,1),:); 
         RP_POS.AvgTitles(i)=RP_POS.All_Titles(AvgInd_cut(i,1));
     end
 end
 
 RP_POS.AvgInd_cut=AvgInd_cut;
 RP_POS.AvgGroups=RP_POS.All_Groups(RP_POS.AvgInd_cut(:,1));

 clear A r i j AvgInd_cut
 
 %% Remove Pooled, QC, Sample 77 
RP_NEG.AvgNorm_wQP=RP_NEG.AvgNorm_cut([14:86 88:end],:);
RP_NEG.AvgTitles_wQP=RP_NEG.AvgTitles([14:86 88:end]);
RP_NEG.AvgGroups_wQP=RP_NEG.AvgGroups([14:86 88:end]);

%% Remove Pooled, QC, Sample 77 
RP_POS.AvgNorm_wQP=RP_POS.AvgNorm_cut([17:89 91:end],:);
RP_POS.AvgTitles_wQP=RP_POS.AvgTitles([17:89 91:end]);
RP_POS.AvgGroups_wQP=RP_POS.AvgGroups([17:89 91:end]);