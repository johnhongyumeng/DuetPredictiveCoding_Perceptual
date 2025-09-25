% This one is to combine all the results table.

clear

tagaddons={'UnpublishednoBH';'PYR_maus1noBH';'PYR_maus2noBH'}; 
% using noBH version to do the BH correction cross mice.

ResultTable_combine=[];
for i_tag =1: length(tagaddons)
    tagaddon=tagaddons{i_tag};
    load(['./Figs/' tagaddon 'TuningCurve.mat'])
    ResultTable_combine=[ResultTable_combine;Results_table];
end


tagaddons={'Unpublished8degree';'PYR_maus18degree';'PYR_maus28degree'};
% using 8degree version to test the tuning curve.

Tuning_response_mat=[];
for i_tag =1: length(tagaddons)
    tagaddon=tagaddons{i_tag};
    load(['./Figs/' tagaddon 'TuningCurve.mat'])
    Tuning_response_mat=[Tuning_response_mat;Tuning_response_post_mat];
end


Results_table=ResultTable_combine;
n_neu=size(Results_table,1);
% BH correction based on all the cells.
for i_phase= 1:2
    for i_col=5:7
        raw_p_values= Results_table(:,i_col,i_phase);
        adj_p_values = mafdr(raw_p_values, 'BHFDR', true); % Apply Benjamini-Hochberg correction
        Results_table(:,i_col,i_phase)=adj_p_values;
    end

    for i_cell= 1:n_neu
        for i_col=5:6   % ctrl and dev column
            if Results_table(i_cell,i_col,i_phase)>=0.05
                Results_table(i_cell,i_col+3,i_phase)=0;  
            elseif Results_table(i_cell,i_col*2-9,i_phase)<Results_table(i_cell,i_col*2-8,i_phase)
                Results_table(i_cell,i_col+3,i_phase)=1;  
            elseif Results_table(i_cell,i_col*2-9,i_phase)>Results_table(i_cell,i_col*2-8,i_phase)
                Results_table(i_cell,i_col+3,i_phase)=-1;
            end
        
        end
        i_col=7;  % ctrl vs dev
        if Results_table(i_cell,i_col,i_phase)>=0.05
            Results_table(i_cell,i_col+3,i_phase)=0;  
        elseif Results_table(i_cell,2,i_phase)<Results_table(i_cell,4,i_phase)
            Results_table(i_cell,i_col+3,i_phase)=1;  
        elseif Results_table(i_cell,2,i_phase)>Results_table(i_cell,4,i_phase)
            Results_table(i_cell,i_col+3,i_phase)=-1;
        end
    end
end




tagaddon='CombinedBH';
save(['./Figs/' tagaddon 'TuningCurve.mat'],'Results_table')


column_names = {'Pre1', 'Ctrl', 'Pre2', 'Dev', ...
                'P1: Ctrl', 'P2: Dev', 'P3: Dev-Ctrl', ...
                'Sig1: Ctrl', 'Sig2: Dev', 'Sig3: Dev-Ctrl'};


% Extract slices from the third dimension
data1 = array2table(Results_table(:,:,1), 'VariableNames', column_names);
data2 = array2table(Results_table(:,:,2), 'VariableNames', column_names);

filename = ['./Figs/' tagaddon 'Results.xlsx']; % Define output file name
% Write the first dataset to "Sheet1"
writetable(data1, filename, 'Sheet', 'odd 1 norm');

% Write the second dataset to "Sheet2"
writetable(data2, filename, 'Sheet', 'odd 2 norm');


column_degs = {'odd 1-45', 'odd 1-22.5', 'odd 1', 'odd 1+22.5', ...
                'odd 1+45', 'odd 2-22.5', 'odd 2', 'odd 2+ 22.5'};

data3 = array2table(Tuning_response_mat, 'VariableNames', column_degs);

writetable(data3, filename, 'Sheet', 'tuning curve');

disp('Excel file successfully created with separate sheets!');

%%% The following part is for plot the tuning curve for each one. 
nPE_collect=[];
pPE_collect=[];

for stim=1:2   % odd 1 dev> control or odd 2 dev > control
    for i_cell=1:n_neu
        if Results_table(i_cell,10,stim)>0
            if stim==1
                temp_tuning_curve=Tuning_response_mat(i_cell,[7,8,1:7]);

            else
                temp_tuning_curve=Tuning_response_mat(i_cell,[3:8,1:3]);

            end
            figure(stim*1000+i_cell)
            plot(temp_tuning_curve)
            if Results_table(i_cell,8,stim)>0
                pPE_collect=[pPE_collect;temp_tuning_curve];
                set(gcf,'units','points','position',[100+i_cell*3,stim*100,400,300])

            else
                nPE_collect=[nPE_collect;temp_tuning_curve];
                set(gcf,'units','points','position',[100+i_cell*3,50+stim*100,400,300])

            end

        end

    end
end
        
figure(5000)

% Define x-tick positions and corresponding labels
x_ticks = [1,5,9];
x_labels = {'-\pi/2', '0', '\pi/2'};
% Define y-tick positions
y_ticks = [-1, 0, 3];

% First subplot: pPE_collect
subplot(1,2,1)
hold on
plot(pPE_collect', 'Color', [0.7 0.7 0.7]) % Plot each row in grey
plot(mean(pPE_collect,1), 'k', 'LineWidth', 2) % Plot mean in black
ylim([-1.5, 4]) % Set y-axis limits
xticks(x_ticks) % Set x-tick positions
yticks(y_ticks) % Set y-ticks at -1, 0, and 3

xticklabels(x_labels) % Set x-tick labels
xlabel('Phase', 'FontSize', 14)
ylabel('z-score', 'FontSize', 14)
% title('average pPE tuning curve', 'FontSize', 16)
set(gca, 'FontSize', 12) % Increase font size for axes
hold off

% Second subplot: nPE_collect
subplot(1,2,2)
hold on
plot(nPE_collect', 'Color', [0.7 0.7 0.7]) % Plot each row in grey
plot(mean(nPE_collect,1), 'k', 'LineWidth', 2) % Plot mean in black
ylim([-1.5, 4]) % Set y-axis limits
xticks(x_ticks) % Set x-tick positions
yticks(y_ticks) % Set y-ticks at -1, 0, and 3

xticklabels(x_labels) % Set x-tick labels
xlabel('Phase', 'FontSize', 14)
ylabel('z-score', 'FontSize', 14)
% title('average nPE tuning curve', 'FontSize', 16)
set(gca, 'FontSize', 12) % Increase font size for axes
hold off


% Preallocate aligned matrix
pPE_aligned = zeros(size(pPE_collect));

% Find the index to align all cells to (e.g., center index)
target_index = ceil(size(pPE_collect, 2) / 2);

for i = 1:size(pPE_collect, 1)
    [~, max_idx] = max(pPE_collect(i, :));
    shift_amt = target_index - max_idx;
    pPE_aligned(i, :) = circshift(pPE_collect(i, :), [0, shift_amt]);
end



% Preallocate aligned matrix
nPE_aligned = zeros(size(nPE_collect));

% Find the index to align all cells to (e.g., center index)
target_index = ceil(size(nPE_collect, 2) / 2);

for i = 1:size(nPE_collect, 1)
    [~, max_idx] = max(nPE_collect(i, :));
    shift_amt = target_index - max_idx;
    nPE_aligned(i, :) = circshift(nPE_collect(i, :), [0, shift_amt]);
end





figure(5001)
% First subplot: pPE_collect
subplot(1,2,1)
hold on
plot(pPE_aligned', 'Color', [0.7 0.7 0.7])
plot(mean(pPE_aligned,1), 'k', 'LineWidth', 2)
ylim([-1.5, 4]) % Set y-axis limits
xticks(x_ticks) % Set x-tick positions
yticks(y_ticks) % Set y-ticks at -1, 0, and 3
xticklabels(x_labels) % Set x-tick labels
xlabel('Phase', 'FontSize', 14)
ylabel('z-score', 'FontSize', 14)
% title('average nPE tuning curve', 'FontSize', 16)
set(gca, 'FontSize', 12) % Increase font size for axes
hold off
% Second subplot: nPE_collect
subplot(1,2,2)
hold on
plot(nPE_aligned', 'Color', [0.7 0.7 0.7]) % Plot each row in grey
plot(mean(nPE_aligned,1), 'k', 'LineWidth', 2) % Plot mean in black
ylim([-1.5, 4]) % Set y-axis limits
xticks(x_ticks) % Set x-tick positions
yticks(y_ticks) % Set y-ticks at -1, 0, and 3

xticklabels(x_labels) % Set x-tick labels
xlabel('Phase', 'FontSize', 14)
ylabel('z-score', 'FontSize', 14)
% title('average nPE tuning curve', 'FontSize', 16)
set(gca, 'FontSize', 12) % Increase font size for axes
hold off



% Load your data
% load('your_data.mat'); % Uncomment this if loading from a .mat file

% Define orientation angles (radians)
orientations = -pi/2 : pi/8 : pi/2-0.1; 

% Number of neurons in each group
nPE_8degree=nPE_collect(:,1:8);
pPE_8degree=pPE_collect(:,1:8);
num_neurons_nPE = size(nPE_8degree, 1);
num_neurons_pPE = size(pPE_8degree, 1);

% Preallocate OSI and CV arrays
OSI_nPE = zeros(num_neurons_nPE, 1);
CV_nPE = zeros(num_neurons_nPE, 1);
ResponsePower_nPE = zeros(num_neurons_nPE, 1);

OSI_pPE = zeros(num_neurons_pPE, 1);
CV_pPE = zeros(num_neurons_pPE, 1);
ResponsePower_pPE = zeros(num_neurons_pPE, 1);

%% Function to calculate OSI and CV
for i = 1:num_neurons_nPE
    responses = nPE_8degree(i, :);
    ResponsePower_nPE(i) = sqrt(sum(responses.^2)); % Compute total response power
    OSI_nPE(i) = compute_osi(responses, orientations);
    % CV_nPE(i) = compute_cv(responses_shifted, orientations);
    CV_nPE(i) = std(responses);

end

for i = 1:num_neurons_pPE
    responses = pPE_8degree(i, :);
    % responses_shifted = responses - min(responses); % Shift responses
    ResponsePower_pPE(i) = sqrt(sum(responses.^2)); % Compute total response power
    OSI_pPE(i) = compute_osi(responses, orientations);
    CV_pPE(i) = std(responses);
end


%% Normalize OSI and CV by Response Power
% OSI_nPE = OSI_nPE .* ResponsePower_nPE;
% CV_nPE = CV_nPE .* ResponsePower_nPE;
% OSI_pPE = OSI_pPE .* ResponsePower_pPE;
% CV_pPE = CV_pPE .* ResponsePower_pPE;


%% Plot Bar Graphs for OSI and CV
figure(6000);

% OSI Plot
subplot(1,3,1);
bar_data = [mean(OSI_nPE), mean(OSI_pPE)];
bar(bar_data, 'FaceColor', 'flat');
hold on;
errorbar(1:2, bar_data, [std(OSI_nPE), std(OSI_pPE)] ./ sqrt([num_neurons_nPE, num_neurons_pPE]), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'nPE', 'pPE'});
ylabel('OSI');
title('Orientation Selectivity Index');
hold off;

% CV Plot
subplot(1,3,2);
bar_data = [mean(CV_nPE), mean(CV_pPE)];
bar(bar_data, 'FaceColor', 'flat');
hold on;
errorbar(1:2, bar_data, [std(CV_nPE), std(CV_pPE)] ./ sqrt([num_neurons_nPE, num_neurons_pPE]), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'nPE', 'pPE'});
ylabel('CV');
title('standard variation');
hold off;

% Mean response
subplot(1,3,3);
bar_data = [mean(ResponsePower_nPE), mean(ResponsePower_pPE)];
bar(bar_data, 'FaceColor', 'flat');
hold on;
errorbar(1:2, bar_data, [std(ResponsePower_nPE), std(ResponsePower_pPE)] ./ sqrt([num_neurons_nPE, num_neurons_pPE]), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'nPE', 'pPE'});
ylabel('Mean power');
title('Mean Response Power');
hold off;


%% Perform t-tests
[~, p_OSI] = ttest2(OSI_nPE, OSI_pPE);
[~, p_CV] = ttest2(CV_nPE, CV_pPE);
[~, p_power] = ttest2(ResponsePower_nPE, ResponsePower_pPE);

% Display Results
fprintf('T-test results:\n');
fprintf('OSI: p = %.4f\n', p_OSI);
fprintf('CV: p = %.4f\n', p_CV);
fprintf('Power: p = %.4f\n', p_power);

figure(6001);

% OSI Plot
bar_data = [mean(OSI_pPE), mean(OSI_nPE)];
b = bar(bar_data, 'FaceColor', 'flat', 'BarWidth', 0.4);  % capture bar object

% Set individual bar colors
b.CData(1,:) = [0.5, 0.5, 0.5];  % grey (R,G,B)
b.CData(2,:) = [0, 0.4470, 0.7410];  % MATLAB's default blue

hold on;
errorbar(1:2, bar_data, [std(OSI_pPE), std(OSI_nPE)] ./ sqrt([num_neurons_pPE, num_neurons_nPE]), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'pPE', 'nPE'});
box off;
hold off;
set(gcf,'units','points','position',[200,200,250,250]);
set(gca,'FontSize',16);
set(gcf,'Color','w');
title(sprintf('Orientation z-score difference (p = %.4f)', p_OSI));
file_tag=['./IndividualFigs/'];
figname= 'TuningCurveOS';
file_name= strcat(file_tag,figname,'.jpg');
saveas(gcf,file_name)        
file_name= strcat(file_tag,figname,'.svg');
saveas(gcf,file_name)        

figure(6002);
% bar_data = [mean(CV_nPE), mean(CV_pPE)];

% std Plot
bar_data = [mean(CV_pPE), mean(CV_nPE)];
b = bar(bar_data, 'FaceColor', 'flat', 'BarWidth', 0.4);  % capture bar object

% Set individual bar colors
b.CData(1,:) = [0.5, 0.5, 0.5];  % grey (R,G,B)
b.CData(2,:) = [0, 0.4470, 0.7410];  % MATLAB's default blue
box off;

hold on;
errorbar(1:2, bar_data, [std(CV_pPE), std(CV_nPE)] ./ sqrt([num_neurons_pPE, num_neurons_nPE]), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'pPE', 'nPE'});
hold off;
set(gcf,'units','points','position',[200,200,250,250]);
set(gca,'FontSize',16);
set(gcf,'Color','w');
title(sprintf('z-score std (p = %.4f)', p_CV));
file_tag=['./IndividualFigs/'];
figname= 'TuningCurveSTD';
file_name= strcat(file_tag,figname,'.jpg');
saveas(gcf,file_name)        
file_name= strcat(file_tag,figname,'.svg');
saveas(gcf,file_name)        

figure(6003);
% bar_data =[mean(ResponsePower_nPE), mean(ResponsePower_pPE)];

% std Plot
bar_data = [mean(ResponsePower_pPE), mean(ResponsePower_nPE)];
b = bar(bar_data, 'FaceColor', 'flat', 'BarWidth', 0.4);  % capture bar object

% Set individual bar colors
b.CData(1,:) = [0.5, 0.5, 0.5];  % grey (R,G,B)
b.CData(2,:) = [0, 0.4470, 0.7410];  % MATLAB's default blue
box off;
hold on;
errorbar(1:2, bar_data, [std(ResponsePower_pPE), std(ResponsePower_nPE)] ./ sqrt([num_neurons_pPE, num_neurons_nPE]), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'pPE', 'nPE'});
hold off;
set(gcf,'units','points','position',[200,200,250,250]);
set(gca,'FontSize',16);
set(gcf,'Color','w');
title(sprintf('z-score ave (p = %.4f)', p_power));
file_tag=['./IndividualFigs/'];
figname= 'TuningCurveAVEZsquare';
file_name= strcat(file_tag,figname,'.jpg');
saveas(gcf,file_name)        
file_name= strcat(file_tag,figname,'.svg');
saveas(gcf,file_name)        
