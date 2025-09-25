% This one is to combine all the results table.

clear

tagaddons={'Unpublished8degree';'PYR_maus18degree';'PYR_maus28degree'};

ResultTable_combine=[];
for i_tag =1: length(tagaddons)
    tagaddon=tagaddons{i_tag};
    load(['./Figs/' tagaddon 'TuningCurve.mat'])
    ResultTable_combine=[ResultTable_combine;Results_table];
end

Tuning_response_mat=[];
for i_tag =1: length(tagaddons)
    tagaddon=tagaddons{i_tag};
    load(['./Figs/' tagaddon 'TuningCurve.mat'])
    Tuning_response_mat=[Tuning_response_mat;Tuning_response_post_mat];
end




Results_table=ResultTable_combine;
tagaddon='Combined';
save(['./Figs/' tagaddon 'TuningCurve.mat'],'Results_table')


column_names = {'Pre1', 'Ctrl', 'Pre2', 'Dev', ...
                'P1: Ctrl', 'P2: Dev', 'P3: Dev-Ctrl', ...
                'Sig1: Ctrl', 'Sig2: Dev', 'Sig3: Dev-Ctrl'};


% Extract slices from the third dimension
data1 = array2table(Results_table(:,:,1), 'VariableNames', column_names);
data2 = array2table(Results_table(:,:,2), 'VariableNames', column_names);

filename = ['./Figs/' tagaddon 'Results.xlsx']; % Define output file name
% Write the first dataset to "Sheet1"
writetable(data1, filename, 'Sheet', '45 degree');

% Write the second dataset to "Sheet2"
writetable(data2, filename, 'Sheet', '135 degree');

disp('Excel file successfully created with separate sheets!');

