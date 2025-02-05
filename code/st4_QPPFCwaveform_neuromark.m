%%clear; clc; close all

%% Setup data input & visualization parameters
dataext='NMv2_HCP_DS'; % extended filename=[data '_' ext];
p2param=['Params_' dataext '.mat']; 
load(['../params/' p2param]);

%%

%% waveform plots

qq = 1;
qpp=QPPs{qq,1};

figure, 

for roi = 1:size(qpp,1)
    plot(squeeze(sum(QPPs{1,1}(roi,:), 1)) / (size(QPPs{1,1}(roi,:), 1))), hold on,

end
xlabel('Time (s)');
ylabel('QPP value');
ylim([-1.5 1.5])
title("All Subdomain Networks")


%%


% Define colors for each network
colors = {
    [0, 1, 1], % CB_rois: Cyan
    [0, 1, 0], % VI_rois: Green
    [0.5, 0.5, 0.5], % PL_rois: Gray
    [1, 0, 1], % SC_rois: Magenta
    [1, 0, 0], % SM_rois: Red
    [0, 0, 1], % HC_rois: Blue
    [1, 0.5, 0] % TN_rois: Orange
};

% ROI indices for each network
CB_rois = 1:13;
VI_rois = 14:25;
PL_rois = 26:36;
SC_rois = 37:54;
SM_rois = 55:68;
HC_rois = 69:90;
TN_rois = 91:105;

% Create new figure
figure;

% Placeholder plots for legend
h = gobjects(7, 1);  % Preallocate graphics object array
h(1) = plot(NaN, NaN, 'Color', colors{1}, 'LineWidth', 2);
hold on;
h(2) = plot(NaN, NaN, 'Color', colors{2}, 'LineWidth', 2);
h(3) = plot(NaN, NaN, 'Color', colors{3}, 'LineWidth', 2);
h(4) = plot(NaN, NaN, 'Color', colors{4}, 'LineWidth', 2);
h(5) = plot(NaN, NaN, 'Color', colors{5}, 'LineWidth', 2);
h(6) = plot(NaN, NaN, 'Color', colors{6}, 'LineWidth', 2);
h(7) = plot(NaN, NaN, 'Color', colors{7}, 'LineWidth', 2);

% Loop through each ROI and plot with network's color
for roi = 1:size(qpp,1)
    % Determine the network the ROI belongs to
    if ismember(roi, CB_rois)
        color = colors{1};
    elseif ismember(roi, VI_rois)
        color = colors{2};
    elseif ismember(roi, PL_rois)
        color = colors{3};
    elseif ismember(roi, SC_rois)
        color = colors{4};
    elseif ismember(roi, SM_rois)
        color = colors{5};
    elseif ismember(roi, HC_rois)
        color = colors{6};
    elseif ismember(roi, TN_rois)
        color = colors{7};
    else
        error('ROI not associated with any network');
    end
    
    % Plot with determined color
    plot(squeeze(sum(QPPs{1,1}(roi,:), 1)) / (size(QPPs{1,1}(roi,:), 1)), 'Color', color, 'LineWidth', 0.5);
    hold on;
end

% Label axes and set legend
xlabel('Time (s)');
ylabel('QPP value');
legend(h, 'CB', 'VI', 'PL', 'SC', 'SM', 'HC', 'TN');
ylim([-1.5 1.5]);



%%

%% Show Major Networks in Brainnetome

CB_rois = 1:13;
VI_rois = 14:25;
PL_rois = 26:36;
SC_rois = 37:54;
SM_rois = 55:68;
HC_rois = 69:90;
TN_rois = 91:105;


CB = squeeze(sum(QPPs{1,1}(CB_rois,:), 1)) / (size(QPPs{1,1}(CB_rois,:), 1));
VI = squeeze(sum(QPPs{1,1}(VI_rois,:), 1)) / (size(QPPs{1,1}(VI_rois,:), 1));
PL = squeeze(sum(QPPs{1,1}(PL_rois,:), 1)) / (size(QPPs{1,1}(PL_rois,:), 1));
SC = squeeze(sum(QPPs{1,1}(SC_rois,:), 1)) / (size(QPPs{1,1}(SC_rois,:), 1));
SM = squeeze(sum(QPPs{1,1}(SM_rois,:), 1)) / (size(QPPs{1,1}(SM_rois,:), 1));
HC = squeeze(sum(QPPs{1,1}(HC_rois,:), 1)) / (size(QPPs{1,1}(HC_rois,:), 1));
TN = squeeze(sum(QPPs{1,1}(TN_rois,:), 1)) / (size(QPPs{1,1}(TN_rois,:), 1));

% string values for graph legends
text1 = "$CB$"; 
text2 = "$VI$";
text3 = "$PL$";
text4 = "$SC$";
text5 = "$SM$";
text6 = "$HC$";
text7 = "$TN$";

waveform = figure;
% Placeholder plots for legend
h(1) = plot(NaN,NaN,'Color',colors{1},'LineWidth',2);
hold on;
h(2) = plot(NaN,NaN,'Color',colors{2},'LineWidth',2);
h(3) = plot(NaN,NaN,'Color',colors{3},'LineWidth',2);
h(4) = plot(NaN,NaN,'Color',colors{4},'LineWidth',2);
h(5) = plot(NaN,NaN,'Color',colors{5},'LineWidth',2);
h(6) = plot(NaN,NaN,'Color',colors{6},'LineWidth',2);
h(7) = plot(NaN,NaN,'Color',colors{7},'LineWidth',2);

plot(CB, 'Color', colors{1}, 'LineWidth', 2.0), hold on,
plot(VI, 'Color', colors{2}, 'LineWidth', 2.0), hold on,
plot(PL, 'Color', colors{3}, 'LineWidth', 2.0), hold on,
plot(SC, 'Color', colors{4}, 'LineWidth', 2.0), hold on,
plot(SM, 'Color', colors{5}, 'LineWidth', 2.0), hold on,
plot(HC, 'Color', colors{6}, 'LineWidth', 2.0), hold on,
plot(TN, 'Color', colors{7}, 'LineWidth', 2.0), hold on,


%h = legend(text1, text2, text3, text4, text5, text6);
%set(h, 'Interpreter', 'latex')
xlabel('Time (s)');
ylabel('QPP value');
%legend(h, 'CB', 'VI', 'PL', 'SC', 'SM', 'HC', 'TN');
ylim([-1 1])

%saveas(waveform, strcat('C:/Users/tjlag/Downloads/', 'hcp_Nsubj_', num2str(subjects), '_TR_', strrep(num2str(TR), '.', ''), '_waveform_1.svg'))


%%

% Combine networks into a matrix
networks = [CB; VI; PL; SC; SM; HC; TN];

% Compute the correlation matrix
corrMatrix = corr(networks');


% Plot the correlation matrix using imagesc
corr_fig = figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Label each network on the axes
networkNames = {'CB', 'VI', 'PL', 'SC', 'SM', 'HC', 'TN'};
set(gca, 'XTick', 1:7, 'XTickLabel', networkNames, 'YTick', 1:7, 'YTickLabel', networkNames, 'TickLabelInterpreter', 'none');
xtickangle(45);

% Display the correlation value in each box
for i = 1:7
    for j = 1:7
        text(j, i, sprintf('%.2f', corrMatrix(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end

title('Correlation Matrix of Domain Networks');



%%


figure,
% Placeholder plots for legend
%h(1) = plot(NaN,NaN,'Color',colors{5},'LineWidth',2); hold on;
%h(2) = plot(NaN,NaN,'Color',colors{6},'LineWidth',2);

plot(HC, 'Color', colors{6}, 'LineWidth', 2.0), hold on,
plot(TN, 'Color', colors{7}, 'LineWidth', 2.0), hold on,


%h = legend(text1, text2, text3, text4, text5, text6);
%set(h, 'Interpreter', 'latex')
xlabel('Time (s)');
ylabel('QPP value');
legend('HC', 'TN');
ylim([-1 1])
title(['HC-TN Correlation: ' num2str(corr(HC', TN'))])


%%

%% Show Major Subdomains in Neuromark

% Define maximally distinct colors for each network
colors = {
    [0, 1, 1],     % CB_CB_rois: Cyan
    [1, 0.5, 0],   % VI_OT_rois: Orange
    [0.5, 0, 0.5], % VI_OC_rois: Purple
    [0.5, 0.5, 0.5], % PL_PL_rois: Gray
    [1, 1, 0],     % SC_ET_rois: Yellow
    [0, 1, 0],     % SC_BG_rois: Green
    [1, 0, 1],     % SM_SM_rois: Magenta
    [0, 0.5, 1],   % HC_IT_rois: Light Blue
    [1, 0.75, 0.8],     % HC_TP_rois: Light Pink
    [0.75, 1, 0.25],     % HC_FR_rois: Lime Green 
    [1, 0, 0], % TN_CE_rois: Red
    [0, 0, 1], % TN_DN_rois: Blue
    [0.25, 0.25, 0.25] % TN_SN_rois: Dark Gray
};


% Define subdomain ROIs
CB_CB_rois = 1:13;
VI_OT_rois = 14:19;
VI_OC_rois = 20:25;
PL_PL_rois = 26:36;
SC_EH_rois = 37:39;
SC_ET_rois = 40:45;
SC_BG_rois = 46:54;
SM_SM_rois = 55:68;
HC_IT_rois = 69:75;
HC_TP_rois = 76:80;
HC_FR_rois = 81:90;
TN_CE_rois = 91:93;
TN_DN_rois = 94:101;
TN_SN_rois = 102:105;

CB_CB = squeeze(sum(QPPs{1,1}(CB_CB_rois,:), 1)) / (size(QPPs{1,1}(CB_CB_rois,:), 1));
VI_OT = squeeze(sum(QPPs{1,1}(VI_OT_rois,:), 1)) / (size(QPPs{1,1}(VI_OT_rois,:), 1));
VI_OC = squeeze(sum(QPPs{1,1}(VI_OC_rois,:), 1)) / (size(QPPs{1,1}(VI_OC_rois,:), 1));
PL_PL = squeeze(sum(QPPs{1,1}(PL_PL_rois,:), 1)) / (size(QPPs{1,1}(PL_PL_rois,:), 1));
SC_ET = squeeze(sum(QPPs{1,1}(SC_ET_rois,:), 1)) / (size(QPPs{1,1}(SC_ET_rois,:), 1));
SC_BG = squeeze(sum(QPPs{1,1}(SC_BG_rois,:), 1)) / (size(QPPs{1,1}(SC_BG_rois,:), 1));
SM_SM = squeeze(sum(QPPs{1,1}(SM_SM_rois,:), 1)) / (size(QPPs{1,1}(SM_SM_rois,:), 1));
HC_IT = squeeze(sum(QPPs{1,1}(HC_IT_rois,:), 1)) / (size(QPPs{1,1}(HC_IT_rois,:), 1));
HC_TP = squeeze(sum(QPPs{1,1}(HC_TP_rois,:), 1)) / (size(QPPs{1,1}(HC_TP_rois,:), 1));
HC_FR = squeeze(sum(QPPs{1,1}(HC_FR_rois,:), 1)) / (size(QPPs{1,1}(HC_FR_rois,:), 1));
TN_CE = squeeze(sum(QPPs{1,1}(TN_CE_rois,:), 1)) / (size(QPPs{1,1}(TN_CE_rois,:), 1));
TN_DN = squeeze(sum(QPPs{1,1}(TN_DN_rois,:), 1)) / (size(QPPs{1,1}(TN_DN_rois,:), 1));
TN_SN = squeeze(sum(QPPs{1,1}(TN_SN_rois,:), 1)) / (size(QPPs{1,1}(TN_SN_rois,:), 1));


% string values for graph legends
text1  = "$CB-CB$"; 
text2  = "$VI-OT$";
text3  = "$VI-OC$";
text4  = "$PL-PL$";
text5  = "$SC-ET$";
text6  = "$SC-BG$";
text7  = "$SM-SM$";
text8  = "$HC-IT$";
text9  = "$HC-TP$";
text10 = "$HC-FR$";
text11 = "$TN-CE$";
text12 = "$TN-DN$";
text13 = "$TN-SN$";

waveform = figure;
% Placeholder plots for legend
h(1) = plot(NaN,NaN,'Color',colors{1},'LineWidth',2);
hold on;
h(2) = plot(NaN,NaN,'Color',colors{2},'LineWidth',2);
h(3) = plot(NaN,NaN,'Color',colors{3},'LineWidth',2);
h(4) = plot(NaN,NaN,'Color',colors{4},'LineWidth',2);
h(5) = plot(NaN,NaN,'Color',colors{5},'LineWidth',2);
h(6) = plot(NaN,NaN,'Color',colors{6},'LineWidth',2);
h(7) = plot(NaN,NaN,'Color',colors{7},'LineWidth',2);
h(8) = plot(NaN,NaN,'Color',colors{8},'LineWidth',2);
h(9) = plot(NaN,NaN,'Color',colors{9},'LineWidth',2);
h(10) = plot(NaN,NaN,'Color',colors{10},'LineWidth',2);
h(11) = plot(NaN,NaN,'Color',colors{11},'LineWidth',2);
h(12) = plot(NaN,NaN,'Color',colors{12},'LineWidth',2);
h(13) = plot(NaN,NaN,'Color',colors{13},'LineWidth',2);


plot(CB_CB, 'Color', colors{1}, 'LineWidth', 2.0), hold on,
plot(VI_OT, 'Color', colors{2}, 'LineWidth', 2.0), hold on,
plot(VI_OC, 'Color', colors{3}, 'LineWidth', 2.0), hold on,
plot(PL_PL, 'Color', colors{4}, 'LineWidth', 2.0), hold on,
plot(SC_ET, 'Color', colors{5}, 'LineWidth', 2.0), hold on,
plot(SC_BG, 'Color', colors{6}, 'LineWidth', 2.0), hold on,
plot(SM_SM, 'Color', colors{7}, 'LineWidth', 2.0), hold on,
plot(HC_IT, 'Color', colors{8}, 'LineWidth', 2.0), hold on,
plot(HC_TP, 'Color', colors{9}, 'LineWidth', 2.0), hold on,
plot(HC_FR, 'Color', colors{10}, 'LineWidth', 2.0), hold on,
plot(TN_CE, 'Color', colors{11}, 'LineWidth', 2.0), hold on,
plot(TN_DN, 'Color', colors{12}, 'LineWidth', 2.0), hold on,
plot(TN_SN, 'Color', colors{13}, 'LineWidth', 2.0), hold on,


%h = legend(text1, text2, text3, text4, text5, text6);
%set(h, 'Interpreter', 'latex')
xlabel('Time (s)');
ylabel('QPP value');
%legend(h, 'CB-CB','VI-OT','VI-OC','PL-PL','SC-ET','SC-BG','SM-SM','HC-IT','HC-TP','HC-FR','TN-CE','TN-DN','TN-SN');
ylim([-1 1])

%saveas(waveform, strcat('C:/Users/tjlag/Downloads/', 'hcp_Nsubj_', num2str(subjects), '_TR_', strrep(num2str(TR), '.', ''), '_waveform_1.svg'))


%%

% Combine networks into a matrix
networks = [CB_CB; VI_OT; VI_OC; PL_PL; SC_ET; SC_BG; SM_SM; HC_IT; HC_TP; HC_FR; TN_CE; TN_DN; TN_SN];

% Compute the correlation matrix
corrMatrix = corr(networks');


% Plot the correlation matrix using imagesc
corr_fig = figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Label each network on the axes
networkNames = {'CB-CB','VI-OT','VI-OC','PL-PL','SC-ET','SC-BG','SM-SM','HC-IT','HC-TP','HC-FR','TN-CE','TN-DN','TN-SN'};
set(gca, 'XTick', 1:13, 'XTickLabel', networkNames, 'YTick', 1:13, 'YTickLabel', networkNames, 'TickLabelInterpreter', 'none');
xtickangle(45);

% Display the correlation value in each box
for i = 1:13
    for j = 1:13
        text(j, i, sprintf('%.2f', corrMatrix(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end

title('Correlation Matrix of Subdomain Networks');


%%




% Example Matrices (Replace with your actual data)
A = corr(QPPs{1,1}', QPPs{1,1}');  % Matrix A (example random data)
B = corr(HCP_NC_long_DS_QPP', HCP_NC_long_DS_QPP');  % Matrix B (example random data)

%A = QPPs{1,1};  % Matrix A (example random data)
%B = HCP_NC_long_DS_QPP;  % Matrix B (example random data)

% Preallocate results matrices
h = zeros(size(A)); % Binary significance matrix
p_values = zeros(size(A)); % P-value matrix

% Perform element-wise paired t-tests
for i = 1:size(A,1)
    for j = 1:size(A,2)
        [h(i,j), p_values(i,j)] = ttest(A(i,j), B(i,j));
    end
end

% Visualization
figure;
imagesc(p_values);  % Plot p-values as a heatmap
colormap(jet);      % Use jet colormap (blue = low p-values, red = high p-values)
colorbar;           % Add color scale
title('Element-wise T-test P-values Heatmap');

% Overlay Significant Values
hold on;
[row, col] = find(h == 1); % Find significant points
scatter(col, row, 100, 'k', 'filled'); % Mark significant points with black dots
hold off;

% Label Axes
xlabel('Column Index');
ylabel('Row Index');






%%





figure,
% Placeholder plots for legend
%h(1) = plot(NaN,NaN,'Color',colors{5},'LineWidth',2); hold on;
%h(2) = plot(NaN,NaN,'Color',colors{6},'LineWidth',2);

plot(TN_CE, 'Color', colors{11}, 'LineWidth', 2.0), hold on,
plot(TN_DN, 'Color', colors{12}, 'LineWidth', 2.0), hold on,


%h = legend(text1, text2, text3, text4, text5, text6);
%set(h, 'Interpreter', 'latex')
xlabel('Time (s)');
ylabel('QPP value');
legend('TN-CE', 'TN-DN');
ylim([-1 1])
title(['Central Network - DMN Correlation: ' num2str(corr(TN_CE', TN_DN'))])


%%

all_rois = [TN_CE_rois TN_DN_rois];

%neuromark_labels = readtable("C:\Users\tjlag\Documents\JupyterNotebooks\keilholz_lab\data\hcp_neuromark\ICN_coordinates_tlagrow_cleaned.csv");
icn_name = neuromark_labels.("individual_label");
ex_index = neuromark_labels.("v2_2_order");


% Create a lookup table for the ICN names
icn_lookup = containers.Map(ex_index, icn_name);

% Extract ICN names for the concatenated ROIs
roi_labels = arrayfun(@(x) icn_lookup(x), all_rois, 'UniformOutput', false);

% Compute the correlation matrix
corr_matrix_a = corr(qpp(all_rois,:)');

figure, 

% Plotting
imagesc(corr_matrix_a);
colorbar; % To display the color scale
title('Subdomanin Correlation Matrix');
xticks(1:size(all_rois,2));
yticks(1:size(all_rois,2));
xticklabels(roi_labels);
yticklabels(roi_labels);
xtickangle(45); % To make the x-axis labels readable
colorbar;
colormap('jet');
caxis([-1, 1]);

% Display the correlation value in each box
for i = 1:size(all_rois,2)
    for j = 1:size(all_rois,2)
        text(j, i, sprintf('%.2f', corr_matrix_a(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end

%%

%neuromark_labels = readtable("C:\Users\tjlag\Documents\JupyterNotebooks\keilholz_lab\data\hcp_neuromark\ICN_coordinates_tlagrow_cleaned.csv");
icn_name = neuromark_labels.("individual_abbrev");
ex_index = neuromark_labels.("v2_2_order");


% Create a lookup table for the ICN names
icn_lookup = containers.Map(ex_index, icn_name);

% Extract ICN names for the concatenated ROIs
roi_labels = arrayfun(@(x) icn_lookup(x), all_rois, 'UniformOutput', false);

all_rois = [1:93 102:105];

% Compute the correlation matrix
corr_matrix = corr([qpp(all_rois,:); TN_DN]');

figure, 

% Plotting
imagesc(corr_matrix);
colorbar; % To display the color scale
title('ROI Correlation Matrix');
xticks(1:size(all_rois,2)+1);
yticks(1:size(all_rois,2)+1);
xticklabels([roi_labels 'TN-DN']);
yticklabels([roi_labels 'TN-DN']);
xtickangle(45); % To make the x-axis labels readable
colorbar;
colormap('jet');
caxis([-1, 1]);

% Display the correlation value in each box
for i = 1:size(all_rois,2)+1
    for j = 1:size(all_rois,2)+1
        text(j, i, sprintf('%.2f', corr_matrix(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end


%%

[~, sort_index] = sort(corr_matrix(end,:));


roi_labels = arrayfun(@(x) icn_lookup(x), all_rois, 'UniformOutput', false);
og_stacked_qpp = qpp(all_rois,:);

% Compute the correlation matrix
corr_matrix_sorted = corr([og_stacked_qpp(sort_index(1:end-1),:); TN_DN]');

figure, 

% Plotting
imagesc(corr_matrix_sorted);
colorbar; % To display the color scale
title('Subnetwork Correlation Matrix');
xticks(1:size(sort_index(1:end-1),2)+1);
yticks(1:size(sort_index(1:end-1),2)+1);
xticklabels([roi_labels(sort_index(1:end-1)) 'TN-DN']);
yticklabels([roi_labels(sort_index(1:end-1)) 'TN-DN']);
xtickangle(45); % To make the x-axis labels readable
colorbar;
colormap('jet');
caxis([-1, 1]);

% Display the correlation value in each box
for i = 1:size(sort_index(1:end-1),2)+1
    for j = 1:size(sort_index(1:end-1),2)+1
        text(j, i, sprintf('%.2f', corr_matrix_sorted(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end

