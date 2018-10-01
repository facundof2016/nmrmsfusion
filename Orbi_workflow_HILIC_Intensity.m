%% Positive mode 
HILIC_POS.Data=HILICPositive;
%% Fix Isotopes and low peak width
j=1;
for i=1:length(HILIC_POS.Data)
    if isnumeric(HILIC_POS.Data{i,4})
        idx(j)=i;
        j=j+1;
    elseif i>3 && ismember(HILIC_POS.Data{i,5}(1),'1')
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans 
HILIC_POS.IsotopeIssues=[HILIC_POS.Data(1:3,:);HILIC_POS.Data(idx,:)];
HILIC_POS.Data(idx,:)=[];
clear idx

j=1;
for i=1:length(HILIC_POS.Data)
    if i>3 && HILIC_POS.Data{i,3}<0.08
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans
HILIC_POS.LowPeakWidth=[HILIC_POS.Data(1:3,:);HILIC_POS.Data(idx,:)];
HILIC_POS.Data(idx,:)=[];
clear idx

%%
HILIC_POS.All_Groups=HILIC_POS.Data(2,5:302);
for i=2:length(HILIC_POS.All_Groups)
    if i<20
        HILIC_POS.All_Groups{i}='Blank';
    elseif i>=20 && i<29
        HILIC_POS.All_Groups{i}='Pooled';
    elseif i>=29 && i<52
        HILIC_POS.All_Groups{i}='QC';
    elseif i==52 || i==53
        HILIC_POS.All_Groups{i}='DDAx';
    elseif i>=54 && i<174
        HILIC_POS.All_Groups{i}='No Recurrence';
    elseif i>=174
        HILIC_POS.All_Groups{i}='Recurrence';
    end
end

HILIC_POS.MASS_RT=cell2mat(HILIC_POS.Data(4:end,[1 2]))';
HILIC_POS.All_Titles=HILIC_POS.Data(3,5:302);
HILIC_POS.All_Titles(find(ismember(HILIC_POS.All_Groups,'DDAx')))=[];
HILIC_POS.All_Titles(find(ismember(HILIC_POS.All_Groups,'Blank')))=[];

HILIC_POS.NormData=cell2mat(HILIC_POS.Data(4:end,5:302))';
HILIC_POS.NormData(find(ismember(HILIC_POS.All_Groups,'DDAx')),:)=[];
HILIC_POS.NormData(find(ismember(HILIC_POS.All_Groups,'Blank')),:)=[];

HILIC_POS.All_Groups(find(ismember(HILIC_POS.All_Groups,'DDAx')))=[];
HILIC_POS.All_Groups(find(ismember(HILIC_POS.All_Groups,'Blank')))=[];

%% Negative mode 
HILIC_NEG.Data=HILICNegative;

%% Fix Isotopes and low peak width
j=1;
for i=1:length(HILIC_NEG.Data)
    if isnumeric(HILIC_NEG.Data{i,4})
        idx(j)=i;
        j=j+1;
    elseif i>3 && ismember(HILIC_NEG.Data{i,5}(1),'1')
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans 
HILIC_NEG.IsotopeIssues=[HILIC_NEG.Data(1:3,:);HILIC_NEG.Data(idx,:)];
HILIC_NEG.Data(idx,:)=[];
clear idx

j=1;
for i=1:length(HILIC_NEG.Data)
    if i>3 &&HILIC_NEG.Data{i,3}<0.08
        idx(j)=i;
        j=j+1;
    end
end
clear i j ans
HILIC_NEG.LowPeakWidth=[HILIC_NEG.Data(1:3,:);HILIC_NEG.Data(idx,:)];
HILIC_NEG.Data(idx,:)=[];
clear idx
%%
HILIC_NEG.All_Groups=HILIC_NEG.Data(2,5:303);
for i=1:length(HILIC_NEG.All_Groups)
    if i<20
        HILIC_NEG.All_Groups{i}='Blank';
    elseif i>=21 && i<29
        HILIC_NEG.All_Groups{i}='Pooled';
    elseif i==52 || i==53
        HILIC_NEG.All_Groups{i}='DDAx';
    elseif i>=29 && i<=51
        HILIC_NEG.All_Groups{i}='QC';
    elseif i>=54 && i<175
        HILIC_NEG.All_Groups{i}='No Recurr';
    elseif i>=175
        HILIC_NEG.All_Groups{i}='Recurr';
    end
end

HILIC_NEG.MASS_RT=cell2mat(HILIC_NEG.Data(4:end,[1 3]))';
HILIC_NEG.All_Titles=HILIC_NEG.Data(3,5:303);
HILIC_NEG.All_Titles(find(ismember(HILIC_NEG.All_Groups,'DDAx')))=[];
HILIC_NEG.All_Titles(find(ismember(HILIC_NEG.All_Groups,'Blank')))=[];

HILIC_NEG.NormData=cell2mat(HILIC_NEG.Data(4:end,5:303))';
HILIC_NEG.NormData(find(ismember(HILIC_NEG.All_Groups,'DDAx')),:)=[];
HILIC_NEG.NormData(find(ismember(HILIC_NEG.All_Groups,'Blank')),:)=[];

HILIC_NEG.All_Groups(find(ismember(HILIC_NEG.All_Groups,'DDAx')))=[];
HILIC_NEG.All_Groups(find(ismember(HILIC_NEG.All_Groups,'Blank')))=[];

%% Remove first few seconds from each of the MS methods
% idx1=find(HILIC_POS.MASS_RT(2,:)>0.75);
idx1=find(HILIC_POS.MASS_RT(2,:)>0.75);

% idx2=find(HILIC_NEG.MASS_RT(2,:)>0.75);
idx2=find(HILIC_NEG.MASS_RT(2,:)>0.75);

HILIC_POS.MASS_RT_cut=HILIC_POS.MASS_RT(:,idx1);
HILIC_NEG.MASS_RT_cut=HILIC_NEG.MASS_RT(:,idx2);

HILIC_POS.NormData_cut=HILIC_POS.NormData(:,idx1);
HILIC_NEG.NormData_cut=HILIC_NEG.NormData(:,idx2);
clear idx1 idx2

idx1=find(HILIC_POS.MASS_RT(2,:)<=0.75);
idx2=find(HILIC_NEG.MASS_RT(2,:)<=0.75);
HILIC_POS.Cut.MASS_RT=HILIC_POS.MASS_RT(:,idx1);
HILIC_NEG.Cut.MASS_RT=HILIC_NEG.MASS_RT(:,idx2);

HILIC_POS.Cut.NormData=HILIC_POS.NormData(:,idx1);
HILIC_NEG.Cut.NormData=HILIC_NEG.NormData(:,idx2);
clear idx2 idx1

%% Average Replicates Negative
%Step 1: Create a matrix with the indeces corresponding to the replicates

AvgInd_cut=zeros(length(HILIC_NEG.All_Titles),3);

for i=1:32 % Pooled and QC
AvgInd_cut(i,1:3)=[i 0 0];
end

HILIC_NEG.All_Titles=strrep(HILIC_NEG.All_Titles, '_1', '');
HILIC_NEG.All_Titles=strrep(HILIC_NEG.All_Titles, '_2', '');
HILIC_NEG.All_Titles=strrep(HILIC_NEG.All_Titles, '_3', '');

%%
j=1; %Samples
for i=33:278
    disp(j)
    if j==13 || j==20 || j==66 || j==69
        AvgInd_cut(i,:)=find(ismember(HILIC_NEG.All_Titles,['S' num2str(j) 'A']));
        if AvgInd_cut(i,1)==(AvgInd_cut(i-1,1))
            AvgInd_cut(i,:)=find(ismember(HILIC_NEG.All_Titles,['S' num2str(j) 'B']));
            j=j+1;
        end
        
    elseif length(find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)])))<3 && ~isempty(find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)])))
        AvgInd_cut(i,1:2)=find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)]));
        j=j+1;
    elseif ~isempty(find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)])))
        AvgInd_cut(i,:)=find(ismember(HILIC_NEG.All_Titles,['S' num2str(j)]));
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
     if AvgInd_cut(i,2)~=0 && AvgInd_cut(i,3)~=0
         HILIC_NEG.AvgNorm_cut(i,:)=mean(HILIC_NEG.NormData_cut(AvgInd_cut(i,:),:),1); 
     elseif AvgInd_cut(i,2)~=0 && AvgInd_cut(i,3)==0
         HILIC_NEG.AvgNorm_cut(i,:)=mean(HILIC_NEG.NormData_cut(AvgInd_cut(i,1:2),:),1); 
     else         
         HILIC_NEG.AvgNorm_cut(i,:)=HILIC_NEG.NormData_cut(AvgInd_cut(i,1),:); 
     end
 end
 
HILIC_NEG.AvgInd_cut=AvgInd_cut;
HILIC_NEG.AvgTitles=HILIC_NEG.All_Titles(HILIC_NEG.AvgInd_cut(:,1));
HILIC_NEG.AvgGroups=HILIC_NEG.All_Groups(HILIC_NEG.AvgInd_cut(:,1));

 clear A r i j AvgInd_cut

   %% Average Replicates Positive
%Step 1: Create a matrix with the indeces corresponding to the replicates

AvgInd_cut=zeros(length(HILIC_POS.All_Titles),3);

for i=1:32 % Pooled and QC
    AvgInd_cut(i,1:3)=[i 0 0];
end

HILIC_POS.All_Titles=strrep(HILIC_POS.All_Titles, '_1', '');
HILIC_POS.All_Titles=strrep(HILIC_POS.All_Titles, '_2', '');
HILIC_POS.All_Titles=strrep(HILIC_POS.All_Titles, '_3', '');

%%
j=1; %Samples
for i=33:277
    disp(j)
    if j==13 || j==20 || j==66 || j==69
         AvgInd_cut(i,:)=find(ismember(HILIC_POS.All_Titles,['S' num2str(j) 'A']));    
        if AvgInd_cut(i,1)==(AvgInd_cut(i-1,1))
            AvgInd_cut(i,:)=find(ismember(HILIC_POS.All_Titles,['S' num2str(j) 'B']));
            j=j+1;
        end
        
    elseif length(find(ismember(HILIC_POS.All_Titles,['S' num2str(j)])))<3 && ~isempty(find(ismember(HILIC_POS.All_Titles,['S' num2str(j)])))
        AvgInd_cut(i,1:2)=find(ismember(HILIC_POS.All_Titles,['S' num2str(j)]));
        j=j+1;
    elseif ~isempty(find(ismember(HILIC_POS.All_Titles,['S' num2str(j)])))
        AvgInd_cut(i,:)=find(ismember(HILIC_POS.All_Titles,['S' num2str(j)]));
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
     if AvgInd_cut(i,2)~=0 && AvgInd_cut(i,3)~=0
         HILIC_POS.AvgNorm_cut(i,:)=mean(HILIC_POS.NormData_cut(AvgInd_cut(i,:),:),1); 
     elseif AvgInd_cut(i,2)~=0 && AvgInd_cut(i,3)==0
         HILIC_POS.AvgNorm_cut(i,:)=mean(HILIC_POS.NormData_cut(AvgInd_cut(i,1:2),:),1); 
     else         
         HILIC_POS.AvgNorm_cut(i,:)=HILIC_POS.NormData_cut(AvgInd_cut(i,1),:); 
     end
 end
 
HILIC_POS.AvgInd_cut=AvgInd_cut;
HILIC_POS.AvgTitles=HILIC_POS.All_Titles(HILIC_POS.AvgInd_cut(:,1));
HILIC_POS.AvgGroups=HILIC_NEG.All_Groups(HILIC_POS.AvgInd_cut(:,1));

 clear A r i j AvgInd_cut
 
 %% Remove Pooled, QC, Sample 77 
HILIC_NEG.AvgNorm_wQP=HILIC_NEG.AvgNorm_cut([33:108 110:end],:);
HILIC_NEG.AvgTitles_wQP=HILIC_NEG.AvgTitles([33:108 110:end]);
HILIC_NEG.AvgGroups_wQP=HILIC_NEG.AvgGroups([33:108 110:end]);

%% Remove Pooled, QC, Sample 77 
HILIC_POS.AvgNorm_wQP=HILIC_POS.AvgNorm_cut([33:108 110:end],:);
HILIC_POS.AvgTitles_wQP=HILIC_POS.AvgTitles([33:108 110:end]);
HILIC_POS.AvgGroups_wQP=HILIC_POS.AvgGroups([33:108 110:end]);