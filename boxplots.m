Classes=HILIC_POS.AvgGroups_wQP;
Class(HILIC_POS.Y==2)={'Remission'};
Class(HILIC_POS.Y==1)={'Recurrence'};

%% boxplot
i=22;
    boxplot(HILIC_POS.AvgNorm_wQP(:,i),Class);
    title([HILIC_POS.ID(i+1,2),'Mass:' num2str(HILIC_POS.MASS_RT(1,i)), 'RT:' num2str(HILIC_POS.MASS_RT(2,i))],'Interpreter','none')
        
    HILIC_POS.Sign.Pvalues(i)

%%
Classes=HILIC_NEG.AvgGroups_wQP;
Class(HILIC_NEG.Y==2)={'Remission'};
Class(HILIC_NEG.Y==1)={'Recurrence'};

%% boxplot
i=27
    boxplot(HILIC_NEG.AvgNorm_wQP(:,i),Class);
    title([HILIC_NEG.ID(i+1,1),'Mass:' num2str(HILIC_NEG.MASS_RT(1,i)), 'RT:' num2str(HILIC_NEG.MASS_RT(2,i))],'Interpreter','none')

    HILIC_NEG.Sign.Pvalues(i)
    
    %% NMR
Classes=NMR.All_Groups;
Class(NMR.Y==2)={'Remission'};
Class(NMR.Y==1)={'Recurrence'};

%% boxplot
i=10
    boxplot(NMR.NormData_sum(:,i),Class);
    title([NMR.ID_sum(i),'Chemical Shift:' num2str(NMR.ppm_sum(i))],'Interpreter','none')

NMR.Sign_sum.Pvalues(i)