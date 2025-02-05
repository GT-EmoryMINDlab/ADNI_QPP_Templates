

%%


% Define the matrices
GrpSn = ADNI_NC_D1; % Replace with your time series matrix
qpp = ADNI_uNC_pre_QPP_rev(:,5:12);    % Replace with your pattern template
trrr = 3;
corr_thresh = 0.2;

% Dimensions
[~, time_length] = size(GrpSn);
[~, template_length] = size(qpp);

% Initialize counters and storage for correlations
count_above_0_2 = 0;
count_below_neg_0_2 = 0;
correlation_values = zeros(1, time_length - template_length + 1);

% Sliding window correlation
for t = 1:(time_length - template_length + 1)
    % Extract the current window from the time series
    current_window = GrpSn(:, t:(t + template_length - 1));
    
    % Compute correlation between the template and the current window
    r = corr(current_window(:), qpp(:));
    correlation_values(t) = r; % Store correlation
    
    % Check the threshold conditions
    if r > corr_thresh
        count_above_0_2 = count_above_0_2 + 1;
    elseif r < -corr_thresh
        count_below_neg_0_2 = count_below_neg_0_2 + 1;
    end
end

% Calculate the total minutes in the time series
total_minutes = (time_length * trrr) / 60;

% Calculate occurrences per minute
occurrences_above_0_2_per_min = count_above_0_2 / total_minutes;
occurrences_below_neg_0_2_per_min = count_below_neg_0_2 / total_minutes;

% Display the results
fprintf('Number of occurrences where correlation > 0.2: %d\n', count_above_0_2);
fprintf('Number of occurrences where correlation < -0.2: %d\n', count_below_neg_0_2);
fprintf('Occurrences per minute where correlation > 0.2: %.2f\n', occurrences_above_0_2_per_min);
fprintf('Occurrences per minute where correlation < -0.2: %.2f\n', occurrences_below_neg_0_2_per_min);

% Plot and save the correlation values
figure;
plot(correlation_values, 'LineWidth', 1.5);
xlabel('Window Index');
ylabel('Correlation');
title('Correlation Values Across Sliding Windows');
grid on;
saveas(gcf, 'correlation_plot.png'); % Save the plot as a PNG

% Plot and save the histogram of correlation values
figure;
histogram(correlation_values, 'BinWidth', 0.05, 'FaceColor', 'blue', 'EdgeColor', 'black');
xlabel('Correlation Value');
ylabel('Frequency');
title('Histogram of Correlation Values');
grid on;
saveas(gcf, 'correlation_histogram.png'); % Save the histogram as a PNG



% Dimensions
[~, time_length] = size(GrpSn);
[~, template_length] = size(qpp);

% Initialize storage for correlations
correlation_values = zeros(1, time_length - template_length + 1);

% Sliding window correlation
for t = 1:(time_length - template_length + 1)
    % Extract the current window from the time series
    current_window = GrpSn(:, t:(t + template_length - 1));
    
    % Compute correlation between the template and the current window
    r = corr(current_window(:), qpp(:));
    correlation_values(t) = r; % Store correlation
end

% Identify indices where correlation exceeds thresholds
above_threshold_indices = find(correlation_values > corr_thresh);
below_threshold_indices = find(correlation_values < -corr_thresh);

% Combine indices and sort
all_threshold_indices = sort([above_threshold_indices, below_threshold_indices]);

% Calculate dwell times between instances
dwell_times = diff(all_threshold_indices) * trrr; % Time gaps in seconds

% Calculate average dwell time
average_dwell_time = mean(dwell_times);

% Display results
fprintf('Average dwell time between instances: %.2f seconds\n', average_dwell_time);
fprintf('Number of dwell time intervals: %d\n', length(dwell_times));

% Plot and save the correlation values
figure;
plot(correlation_values, 'LineWidth', 1.5);
xlabel('Window Index');
ylabel('Correlation');
title('Correlation Values Across Sliding Windows');
grid on;
saveas(gcf, 'correlation_plot_with_dwell.png'); % Save the plot as a PNG

% Plot histogram of dwell times
figure;
histogram(dwell_times, 'BinWidth', 1, 'FaceColor', 'blue', 'EdgeColor', 'black');
xlabel('Dwell Time Between Instances (seconds)');
ylabel('Frequency');
title('Histogram of Dwell Times Between Instances');
set(gca, 'YScale', 'log'); % Set y-axis to log scale
grid on;
saveas(gcf, 'dwell_time_histogram.png'); % Save the histogram as a PNG


%%

% Define network ROIs
CB_rois = 1:13;
VI_rois = 14:25;
PL_rois = 26:36;
SC_rois = 36:54;
SM_rois = 55:68;
HC_rois = 69:90;
TN_rois = 91:105;

% define QPP templates that you wish to compare
QPP_templates = {HCP_NC_long_QPP, ADNI_DAT_US_QPP, ADNI_MCI_US_QPP, ADNI_NC_US_QPP};

% Predefine storage for all networks across templates
allNetworks = [];
allLabels = [];
networkNames = {'CB', 'VI', 'PL', 'SC', 'SM', 'HC', 'TN'};

% Loop through each QPP template
for t = 1:length(QPP_templates)
    % Extract the current template
    currentQPP = QPP_templates{t};

    % Perform domain-wise summation and normalization
    CB = squeeze(sum(currentQPP(CB_rois, :), 1)) / size(currentQPP(CB_rois, :), 1);
    VI = squeeze(sum(currentQPP(VI_rois, :), 1)) / size(currentQPP(VI_rois, :), 1);
    PL = squeeze(sum(currentQPP(PL_rois, :), 1)) / size(currentQPP(PL_rois, :), 1);
    SC = squeeze(sum(currentQPP(SC_rois, :), 1)) / size(currentQPP(SC_rois, :), 1);
    SM = squeeze(sum(currentQPP(SM_rois, :), 1)) / size(currentQPP(SM_rois, :), 1);
    HC = squeeze(sum(currentQPP(HC_rois, :), 1)) / size(currentQPP(HC_rois, :), 1);
    TN = squeeze(sum(currentQPP(TN_rois, :), 1)) / size(currentQPP(TN_rois, :), 1);

    % Combine networks into a matrix
    networks = [CB; VI; PL; SC; SM; HC; TN];

    % Append to the overall networks matrix
    allNetworks = [allNetworks; networks];

    % Append corresponding labels
    allLabels = [allLabels; strcat(networkNames, sprintf('_T%d', t))'];
end

% Compute the correlation matrix for all networks
corrMatrix = corr(allNetworks');

% Plot the correlation matrix using imagesc
figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Label each network on the axes
set(gca, 'XTick', 1:size(allNetworks, 1), 'XTickLabel', allLabels, 'YTick', 1:size(allNetworks, 1), 'YTickLabel', allLabels, 'TickLabelInterpreter', 'none');
xtickangle(45);

% Display the correlation value in each box
for i = 1:size(allNetworks, 1)
    for j = 1:size(allNetworks, 1)
        text(j, i, sprintf('%.2f', corrMatrix(i, j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end

% Add title
title('Correlation Matrix of Domain Networks Across Templates');


%%

% Define network ROIs
CB_rois = 1:13;
VI_rois = 14:25;
PL_rois = 26:36;
SC_rois = 37:54;
SM_rois = 55:68;
HC_rois = 69:90;
TN_rois = 91:105;

% QPP templates
%QPP_templates = {ADNI_ALL_DATA_CONCAT_QPP, ADNI_uNC_pre_QPP, ADNI_NC_QPP, ADNI_uNC_post_QPP, ADNI_pMCI_pre_QPP, ADNI_MCI_QPP, ADNI_pMCI_post_QPP, ADNI_DAT_QPP};
%QPP_templates = {ADNI_NC_QPP, ADNI_MCI_QPP, ADNI_DAT_QPP, ADNI_uNC_pre_QPP_rev, ADNI_uNC_post_QPP, ADNI_pMCI_pre_QPP, ADNI_pMCI_post_QPP, };
QPP_templates = {ADNI_NC_QPP};

% Concatenate all networks from each template
allNetworks = [];

for t = 1:length(QPP_templates)
    % Extract the current template
    currentQPP = QPP_templates{t}';

    % Append the entire network data for the template
    allNetworks = [allNetworks, currentQPP];
end

% Compute the correlation matrix for all networks
corrMatrix = corr(allNetworks);

% Plot the correlation matrix using imagesc
figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Add lines between templates
hold on;
templateSize = size(QPP_templates{1}, 1);
for t = 1:length(QPP_templates) - 1
    linePos = t * templateSize;
    xline(linePos, 'k-', 'LineWidth', 5);
    yline(linePos, 'k-', 'LineWidth', 5);
end
axis square;
hold off;

% Add title
title('Correlation Matrix of All 105 Networks Across Templates');


%%


% Example Matrices (Replace with your actual data)
A = corr(ADNI_NC_QPP', ADNI_NC_QPP');  % Matrix A (example random data)
B = corr(ADNI_DAT_QPP', ADNI_DAT_QPP');  % Matrix B (example random data)

% Assuming matrices A and B of the same size (NxN)
diff_mat = A - B;
p_values = arrayfun(@(i) signrank(A(i), B(i)), 1:numel(A));
p_values = reshape(p_values, size(A));



figure;
heatmap(diff_mat);
title('Element-wise Difference Heatmap');


%%

A = corr(ADNI_uNC_pre_QPP(:,5:12)', ADNI_uNC_pre_QPP(:,5:12)');
B = corr(ADNI_pMCI_post_QPP(:,5:12)', ADNI_pMCI_post_QPP(:,5:12)');
A_minus_B = A-B;



figure,

subplot(1,3,1), 
imagesc(A), 
axis square;
axis off;

subplot(1,3,2), 
imagesc(B),
axis square;
axis off;

subplot(1,3,3),

imagesc(A_minus_B), 


axis square;
axis off;


colormap('jet')


% Assuming you have a vector of variances, one per network
network_variances = var(A_minus_B, 0, 2); % Row-wise variance of A - B
expected_variance = mean(network_variances); % Null hypothesis: Variances follow an expected mean

% Chi-square test per network
df = 104; % Degrees of freedom (N-1 where N=105 samples)
chi_sq_values = (network_variances / expected_variance) * df; % Compute chi-square statistic
p_values = 1 - chi2cdf(chi_sq_values, df); % Right-tailed test

% Find significant networks at p < 0.05
significant_networks = find(p_values < 0.05);

% Display results
disp('Significant Networks:');
disp(significant_networks);


% Define domains
CB_rois = 1:13;
VI_rois = 14:25;
PL_rois = 26:36;
SC_rois = 37:54;
SM_rois = 55:68;
HC_rois = 69:90;
TN_rois = 91:105;

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

% Group all domains and subdomains
roi_domains = {CB_rois, VI_rois, PL_rois, SC_rois, SM_rois, HC_rois, TN_rois};
roi_subdomains = {CB_CB_rois, VI_OT_rois, VI_OC_rois, PL_PL_rois, SC_EH_rois, SC_ET_rois, SC_BG_rois, ...
                  SM_SM_rois, HC_IT_rois, HC_TP_rois, HC_FR_rois, TN_CE_rois, TN_DN_rois, TN_SN_rois};

% Define domain boundaries (used for bar coloring)
domain_boundaries = [13, 25, 36, 54, 68, 90, 105];

% Define subdomain boundaries (used for vertical lines)
subdomain_boundaries = [13, 19, 25, 36, 39, 45, 54, 68, 75, 80, 90, 93, 101, 105];

% Define custom colors for each domain (normalized RGB values)
custom_colors = [
    0, 255, 255;  % Light Blue (CB)
    0, 255, 0;  % Purple (VI)
    152, 152,  152;  % Yellow (PL)
    255, 0, 255;  % Blue-Green (SC)
    255, 0, 0;  % Green-Purple (SM)
    0, 0, 255;  % Pink-Purple (HC)
    255, 128, 0   % Pink-Yellow (TN)
] / 255; % Normalize to [0,1] range for MATLAB

% Create bar graph
figure;
b = bar(network_variances, 'FaceColor', 'flat'); % Enable individual coloring
hold on;

% Assign colors to bars based on their domain
for i = 1:length(roi_domains)
    idx = roi_domains{i}; % Get indices of current domain
    b.CData(idx, :) = repmat(custom_colors(i, :), length(idx), 1); % Assign color
end

% Highlight significant networks
scatter(significant_networks, network_variances(significant_networks), 75, 'r', 'filled'); 

% Add vertical black lines at subdomain boundaries
for i = 1:length(subdomain_boundaries)
    xline(subdomain_boundaries(i) + 0.5, 'k', 'LineWidth', 1.5); % Align between bars
end

% Formatting
axis square;
xlim([0, 105]);
ylim([0, 1.5]);

hold off;





%%



QPP_templates = {ADNI_NC_QPP, ADNI_MCI_QPP, ADNI_DAT_QPP, ADNI_uNC_pre_QPP_rev, ADNI_uNC_post_QPP, ADNI_pMCI_pre_QPP, ADNI_pMCI_post_QPP};

% Concatenate all networks from each template
allNetworks = [];

for t = 1:length(QPP_templates)
    % Extract the current template
    currentQPP = QPP_templates{t}';

    % Append the entire network data for the template
    allNetworks = [allNetworks, currentQPP];
end

% Compute the correlation matrix for all networks
corrMatrix = corr(allNetworks);

% Modify the matrix: lower triangular = correlation, upper triangular = absolute correlation
absCorrMatrix = abs(corrMatrix); % Absolute correlation matrix
for i = 1:size(corrMatrix, 1)
    for j = i+1:size(corrMatrix, 2)
        corrMatrix(i, j) = absCorrMatrix(i, j); % Set upper triangular to absolute values
        corrMatrix(j, i) = corrMatrix(j, i); % Keep lower triangular as raw correlation
    end
end

% Plot the modified correlation matrix using imagesc
figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Add lines between templates
hold on;
templateSize = size(QPP_templates{1}, 1);
for t = 1:length(QPP_templates) - 1
    linePos = t * templateSize;
    xline(linePos, 'k-', 'LineWidth', 4);
    yline(linePos, 'k-', 'LineWidth', 4);
end
hold off;
axis square;
axis off;
% Add title
%title('Correlation Matrix with Lower = Correlation, Upper = Absolute Correlation');


%%

% Define network ROIs
CB_rois = 1:13;
VI_rois = 14:25;
PL_rois = 26:36;
SC_rois = 37:54;
SM_rois = 55:68;
HC_rois = 69:90;
TN_rois = 91:105;


%QPP_templates = {ADNI_NC_QPP, ADNI_MCI_QPP, ADNI_DAT_QPP, ADNI_uNC_pre_QPP_rev, ADNI_uNC_post_QPP, ADNI_pMCI_pre_QPP, ADNI_pMCI_post_QPP};
QPP_templates = {ADNI_NC_QPP};

% Concatenate all networks from each template
allNetworks = [];

for t = 1:length(QPP_templates)
    % Extract the current template
    currentQPP = QPP_templates{t}';

    % Append the entire network data for the template
    allNetworks = [allNetworks, currentQPP];
end

% Compute the correlation matrix for all networks
corrMatrix = corr(allNetworks);

% Modify the matrix: lower triangular = correlation, upper triangular = absolute correlation
absCorrMatrix = abs(corrMatrix); % Absolute correlation matrix



for i = 1:size(corrMatrix, 1)
    for j = i+1:size(corrMatrix, 2)
        corrMatrix(i, j) = absCorrMatrix(i, j); % Set upper triangular to absolute values
        corrMatrix(j, i) = corrMatrix(j, i); % Keep lower triangular as raw correlation
    end
end



% Plot the modified correlation matrix using imagesc
figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Add lines between templates
hold on;
templateSize = size(QPP_templates{1}, 1);
for t = 1:length(QPP_templates) - 1
    linePos = t * templateSize;
    xline(linePos, 'k-', 'LineWidth', 4);
    yline(linePos, 'k-', 'LineWidth', 4);
end
hold off;

hold off;



%%

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

% Define all subdomains
subdomainRanges = {CB_CB_rois, VI_OT_rois, VI_OC_rois, PL_PL_rois, ...
    SC_EH_rois, SC_ET_rois, SC_BG_rois, SM_SM_rois, HC_IT_rois, ...
     HC_TP_rois, HC_FR_rois, TN_CE_rois, TN_DN_rois, TN_SN_rois};

% Define network ROIs
CB_rois = 1:13;
VI_rois = 14:25;
PL_rois = 26:36;
SC_rois = 37:54;
SM_rois = 55:68;
HC_rois = 69:90;
TN_rois = 91:105;

% Define abbreviations for each network
networkNames = {'CB', 'VI', 'PL', 'SC', 'SM', 'HC', 'TN'};
networkRanges = {CB_rois, VI_rois, PL_rois, SC_rois, SM_rois, HC_rois, TN_rois};

% Concatenate all networks from each template
QPP_templates = {ADNI_pMCI_post_QPP(:,5:12)};
allNetworks = [];

for t = 1:length(QPP_templates)
    % Extract the current template
    currentQPP = QPP_templates{t}';
    % Append the entire network data for the template
    allNetworks = [allNetworks, currentQPP];
end

% Compute the correlation matrix for all networks
corrMatrix = corr(allNetworks);

% Modify the matrix: lower triangular = correlation, upper triangular = absolute correlation
absCorrMatrix = abs(corrMatrix); % Absolute correlation matrix
for i = 1:size(corrMatrix, 1)
    for j = i+1:size(corrMatrix, 2)
        corrMatrix(i, j) = corrMatrix(i, j); % Set upper triangular to absolute values
        corrMatrix(j, i) = corrMatrix(j, i); % Keep lower triangular as raw correlation
    end
end

% Plot the modified correlation matrix using imagesc
figure;
imagesc(A_minus_B);
colorbar;
colormap('jet');
caxis([-1, 1]);
hold on;

% Add labels for networks
for k = 1:length(networkRanges)
    range = networkRanges{k};
    middleIndex = mean(range); % Compute the middle index of the network range
    % Add labels at the middle of the range for both x and y axes
    text(middleIndex, -5, networkNames{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold'); % Bottom x-axis
    text(-5, middleIndex, networkNames{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90); % Left y-axis
end

% Add lines at the borders of subdomains
for k = 1:length(subdomainRanges)
    range = subdomainRanges{k};
    borderStart = min(range) - 0.5; % Subtract 0.5 for proper alignment with grid
    borderEnd = max(range) + 0.5; % Add 0.5 for proper alignment with grid

    % Horizontal and vertical lines
    line([borderStart, borderStart], [0.5, size(corrMatrix, 1) + 0.5], 'Color', 'k', 'LineWidth', 2.5); % Left vertical line
    line([borderEnd, borderEnd], [0.5, size(corrMatrix, 1) + 0.5], 'Color', 'k', 'LineWidth', 2.5); % Right vertical line
    line([0.5, size(corrMatrix, 2) + 0.5], [borderStart, borderStart], 'Color', 'k', 'LineWidth', 2.5); % Top horizontal line
    line([0.5, size(corrMatrix, 2) + 0.5], [borderEnd, borderEnd], 'Color', 'k', 'LineWidth', 2.5); % Bottom horizontal line
end

axis square;
axis off;
hold off;



%%

% Define the QPP templates
QPP_templates = {ADNI_NC_QPP, ADNI_MCI_QPP, ADNI_DAT_QPP, ...
    ADNI_uNC_pre_QPP_rev, ADNI_uNC_post_QPP, ADNI_pMCI_pre_QPP, ADNI_pMCI_post_QPP};

% Number of templates
numTemplates = length(QPP_templates);

% Initialize the correlation matrix
corrMatrix = zeros(numTemplates, numTemplates);


% Compute correlations
for i = 1:numTemplates
    for j = 1:numTemplates
        % Calculate the correlation between templates
        corrValue = mean(mean(abs(corr(QPP_templates{i}', QPP_templates{j}'))));
        corrMatrix(i, j) = corrValue;
    end
end

% Display the correlation matrix
disp('Correlation Matrix:');
disp(corrMatrix);

% Define prognosis labels
prognosisLabels = {'sNC', 'sMCI', 'sDAT', 'uNC-pre', 'uNC-post', 'pMCI-pre', 'pMCI-post'};

% Calculate the self-correlations (diagonal of the matrix)
selfCorr = diag(corrMatrix);

% Create a matrix where self-correlation values are subtracted from the rest
selfCorrMatrix = corrMatrix; % Initialize with the original
for i = 1:size(corrMatrix, 1)
    selfCorrMatrix(i, :) = selfCorr(i) - corrMatrix(i, :);
end

% Plotting the results
figure;
%{
% First subplot: Regular correlation matrix
subplot(1, 2, 1);
imagesc(corrMatrix);
colormap(jet); % Set color scheme to jet
colorbar; % Add color bar to show the scale
title('Correlation Matrix');
xlabel('Templates');
ylabel('Templates');
axis square; % Make it square
set(gca, 'XTick', 1:length(prognosisLabels), 'XTickLabel', prognosisLabels, ...
         'YTick', 1:length(prognosisLabels), 'YTickLabel', prognosisLabels); % Set labels
caxis([0, 1]); % Set color scale from 0 to 1


%}


% Second subplot: Self-correlation minus present value
%subplot(1, 2, 2);
imagesc(selfCorrMatrix);
colormap("autumn"); % Set color scheme to jet
colorbar; % Add color bar to show the scale
%title('Self-Corr Minus Present Value');
%xlabel('Templates');
%ylabel('Templates');
%axis square; % Make it square
set(gca, 'XTick', 1:length(prognosisLabels), 'XTickLabel', prognosisLabels, ...
         'YTick', 1:length(prognosisLabels), 'YTickLabel', prognosisLabels); % Set labels
caxis([0, .3]); % Set color scale from 0 to 1



%%

% Define the QPP templates
QPP_templates = {ADNI_NC_QPP, ADNI_MCI_QPP, ADNI_DAT_QPP, ...
    ADNI_uNC_pre_QPP_rev, ADNI_uNC_post_QPP, ADNI_pMCI_pre_QPP, ADNI_pMCI_post_QPP};

% Number of templates
numTemplates = length(QPP_templates);

% Initialize the 3D correlation matrices
fullCorrMatrices = cell(numTemplates, numTemplates);
selfCorrMatrices = cell(numTemplates, 1);

% Compute the full correlation matrices
for i = 1:numTemplates
    % Self-correlation for each template
    selfCorrMatrices{i} = corr(QPP_templates{i}(:,5:12)', QPP_templates{i}(:,5:12)');
    for j = 1:numTemplates
        % Correlation between templates
        fullCorrMatrices{i, j} = corr(QPP_templates{i}(:,5:12)', QPP_templates{j}(:,5:12)');
    end
end

% Initialize the difference matrices
diffMatrices = cell(numTemplates, numTemplates);

% Compute the difference from the self-correlation matrix
for i = 1:numTemplates
    for j = 1:numTemplates
        diffMatrices{i, j} = selfCorrMatrices{i} - fullCorrMatrices{i, j};
    end
end

% Define prognosis labels
prognosisLabels = {'sNC', 'sMCI', 'sDAT', 'uNC-pre', 'uNC-post', 'pMCI-pre', 'pMCI-post'};

% Visualize the results
figure;

% Plot example: full correlation matrix (sNC with sMCI)
subplot(1, 2, 1);
imagesc(fullCorrMatrices{1, 2}); % Change indices to visualize specific matrices
colormap(jet);
colorbar;
title('Full Correlation Matrix (sNC vs sMCI)');
xlabel('Timepoints/Features');
ylabel('Timepoints/Features');
axis square;

% Plot example: difference matrix (sNC self-correlation minus sNC vs sMCI)
subplot(1, 2, 2);
imagesc(diffMatrices{1, 2}); % Change indices to visualize specific matrices
colormap(jet);
colorbar;
title('Difference Matrix (Self Corr - sNC vs sMCI)');
xlabel('Timepoints/Features');
ylabel('Timepoints/Features');
axis square;


%%

% Define the size of the individual matrices
matrixSize = size(diffMatrices{1, 1}, 1); % Assuming all are 105x105

% Create a large composite matrix to hold the 7x7 grid of 105x105 matrices
compositeMatrix = zeros(matrixSize * numTemplates, matrixSize * numTemplates);

% Fill the composite matrix with the 7x7 diffMatrices
for i = 1:numTemplates
    for j = 1:numTemplates
        % Calculate the row and column ranges for the current matrix
        rowStart = (i - 1) * matrixSize + 1;
        rowEnd = i * matrixSize;
        colStart = (j - 1) * matrixSize + 1;
        colEnd = j * matrixSize;

        % Place the current difference matrix in the composite matrix
        compositeMatrix(rowStart:rowEnd, colStart:colEnd) = diffMatrices{i, j};
    end
end

% Visualization
figure;
imagesc(compositeMatrix);
colormap(jet);
colorbar;
%title('Composite Plot of 7x7 Difference Matrices');
%xlabel('Templates');
%ylabel('Templates');
%axis square;

% Add gridlines to separate the 105x105 matrices
hold on;
for k = 1:numTemplates
    % Vertical lines
    xline(k * matrixSize + 0.5, 'k', 'LineWidth', 2);
    % Horizontal lines
    yline(k * matrixSize + 0.5, 'k', 'LineWidth', 2);
end

axis square;
axis off
hold off;



%%

% Example matrices (replace with your actual matrices)
matrix1 = abs(fullCorrMatrices{1,4}); % Matrix 1
matrix2 = abs(fullCorrMatrices{1,5}); % Matrix 2
matrix3 = abs(fullCorrMatrices{1,6}); % Matrix 2
matrix4 = abs(fullCorrMatrices{1,7}); % Matrix 2


% Generate example data
group1 = mean(matrix1,2); % Replace with your actual group data
group2 = mean(matrix2,2); % Replace with your actual group data
group3 = mean(matrix3,2); % Replace with your actual group data
group4 = mean(matrix4,2); % Replace with your actual group data

% Combine data and create group labels
data = [group1; group2; group3; group4];
groups = [ones(105, 1); 2 * ones(105, 1); 3 * ones(105, 1); 4 * ones(105, 1)];

% Perform Kruskal-Wallis test
[p, tbl, stats] = kruskalwallis(data, groups, 'off'); % 'off' suppresses the plot
disp(['Kruskal-Wallis p-value: ', num2str(p)]);



% Flatten the matrices into vectors for pairwise testing
values1 = matrix1(:);
values2 = matrix2(:);
values3 = matrix3(:);
values4 = matrix4(:);

% Combine the data and create group labels
data = [values1; values2; values3; values4];
groups = [ones(length(values1), 1); 2 * ones(length(values2), 1); 3 * ones(length(values3), 1); 4 * ones(length(values4), 1)];

% Perform Kruskal-Wallis Test
[p, tbl, stats] = kruskalwallis(data, groups, 'off');
disp(['Kruskal-Wallis Test p-value: ', num2str(p)]);

% Boxplot for visual comparison
boxplot(data, groups, 'Labels', {'Matrix 1', 'Matrix 2', 'Matrix 3', 'Matrix 4'});
ylabel('Values');
title('Comparison of Matrix Distributions');


%%

% Example matrices (replace with your actual matrices)
matrix1 = fullCorrMatrices{3,1}; % Matrix 1
matrix2 = fullCorrMatrices{3,2}; % Matrix 2
matrix3 = fullCorrMatrices{3,3}; % Matrix 3

% Generate example data
group1 = mean(matrix1,2); % Replace with your actual group data
group2 = mean(matrix2,2); % Replace with your actual group data
group3 = mean(matrix3,2); % Replace with your actual group data

% Combine data and create group labels
data = [group1; group2; group3];
groups = [ones(length(group1), 1); 2 * ones(length(group2), 1); ...
          3 * ones(length(group3), 1)];



% Perform Kruskal-Wallis test
[p, tbl, stats] = kruskalwallis(data, groups, 'off'); % 'off' suppresses the plot
disp(['Kruskal-Wallis p-value: ', num2str(p)]);

% Boxplot for visual comparison
figure;
boxplot(data, groups, 'Labels', {'sNC', 'sMCI', 'sDAT'});
ylabel('Values');
title('Comparison of Matrix Distributions');

% Post-hoc analysis: Dunn test
disp('--- Post-Hoc Analysis: Dunn Test ---');
comparisons = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');
disp('Pairwise comparisons using Dunn Test:');
disp(array2table(comparisons, 'VariableNames', ...
    {'Group 1', 'Group 2', 'Lower Bound', 'Difference', 'Upper Bound', 'p-value'}));

% Post-hoc analysis: Pairwise Mann-Whitney U Tests
disp('--- Post-Hoc Analysis: Pairwise Mann-Whitney U Tests ---');
pairs = nchoosek(1:3, 2); % Generate all pairwise combinations of groups
for i = 1:size(pairs, 1)
    groupA = data(groups == pairs(i, 1)); % Group 1
    groupB = data(groups == pairs(i, 2)); % Group 2
    p_mwu = ranksum(groupA, groupB); % Mann-Whitney U test
    fprintf('Comparison between Group %d and Group %d: p-value = %.4f\n', ...
            pairs(i, 1), pairs(i, 2), p_mwu);
end



%%

% Example matrices (replace with your actual matrices)
matrix1 = abs(fullCorrMatrices{1,4}); % Matrix 1
matrix2 = abs(fullCorrMatrices{1,5}); % Matrix 2
matrix3 = abs(fullCorrMatrices{1,6}); % Matrix 3
matrix4 = abs(fullCorrMatrices{1,7}); % Matrix 4

% Generate example data
group1 = mean(matrix1,2); % Replace with your actual group data
group2 = mean(matrix2,2); % Replace with your actual group data
group3 = mean(matrix3,2); % Replace with your actual group data
group4 = mean(matrix4,2); % Replace with your actual group data

% Combine data and create group labels
data = [group1; group2; group3; group4];
groups = [ones(length(group1), 1); 2 * ones(length(group2), 1); ...
          3 * ones(length(group3), 1); 4 * ones(length(group4), 1)];



% Perform Kruskal-Wallis test
[p, tbl, stats] = kruskalwallis(data, groups, 'off'); % 'off' suppresses the plot
disp(['Kruskal-Wallis p-value: ', num2str(p)]);

% Boxplot for visual comparison
figure;
boxplot(data, groups, 'Labels', {'uNC-pre', 'uNC-post', 'pMCI-pre', 'pMCI-post'});
ylabel('Values');
title('Comparison of Matrix Distributions');

% Post-hoc analysis: Dunn test
disp('--- Post-Hoc Analysis: Dunn Test ---');
comparisons = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');
disp('Pairwise comparisons using Dunn Test:');
disp(array2table(comparisons, 'VariableNames', ...
    {'Group 1', 'Group 2', 'Lower Bound', 'Difference', 'Upper Bound', 'p-value'}));

% Post-hoc analysis: Pairwise Mann-Whitney U Tests
disp('--- Post-Hoc Analysis: Pairwise Mann-Whitney U Tests ---');
pairs = nchoosek(1:4, 2); % Generate all pairwise combinations of groups
for i = 1:size(pairs, 1)
    groupA = data(groups == pairs(i, 1)); % Group 1
    groupB = data(groups == pairs(i, 2)); % Group 2
    p_mwu = ranksum(groupA, groupB); % Mann-Whitney U test
    fprintf('Comparison between Group %d and Group %d: p-value = %.4f\n', ...
            pairs(i, 1), pairs(i, 2), p_mwu);
end


%%


% Define categories
categories = {'sNC', 'sMCI', 'sDAT', 'uNC-PRE', 'uNC-POST', 'pMCI-PRE', 'pMCI-POST'};

% Define mean values
mean_values = [
    0.000, 0.033, 0.119, 0.149, 0.203, 0.141, 0.093;
    0.067, 0.000, 0.080, 0.118, 0.191, 0.122, 0.127;
    0.116, 0.042, 0.000, 0.047, 0.094, 0.068, 0.197;
    0.164, 0.099, 0.065, 0.000, 0.074, 0.124, 0.228;
    0.185, 0.138, 0.079, 0.040, 0.000, 0.105, 0.237;
    0.138, 0.084, 0.067, 0.106, 0.120, 0.000, 0.237;
    0.066, 0.066, 0.174, 0.187, 0.228, 0.214, 0.000
];

% Define standard deviations
std_values = [
    0.000, 0.064, 0.076, 0.113, 0.090, 0.088, 0.112;
    0.069, 0.000, 0.127, 0.144, 0.109, 0.147, 0.095;
    0.090, 0.131, 0.000, 0.139, 0.119, 0.097, 0.083;
    0.114, 0.120, 0.122, 0.000, 0.100, 0.093, 0.090;
    0.077, 0.096, 0.109, 0.091, 0.000, 0.111, 0.064;
    0.081, 0.103, 0.085, 0.100, 0.115, 0.000, 0.064;
    0.108, 0.112, 0.082, 0.060, 0.069, 0.057, 0.000
];

% Define significance matrix (1 for significant, 0 for not significant)
significant_values = [
    0, 0, 1, 1, 1, 1, 1;
    0, 0, 1, 1, 1, 1, 1;
    1, 1, 0, 1, 1, 1, 0;
    1, 1, 1, 0, 0, 1, 1;
    1, 1, 1, 0, 0, 1, 1;
    1, 1, 1, 1, 1, 0, 1;
    1, 1, 0, 1, 1, 1, 0
];

% Define custom colors (normalized RGB values)
custom_colors = [
    172, 249, 255;  % Light Blue (sNC)
    218, 177, 255;  % Purple (sMCI)
    255, 255,  32;  % Yellow (sDAT)
    199, 247, 222;  % Blue-Green (uNC-pre)
    214, 212, 233;  % Green-Purple (uNC-post)
    248, 196, 242;  % Pink-Purple (pMCI-pre)
    254, 220, 124   % Pink-Yellow (pMCI-post)
] / 255; % Normalize to [0,1] range for MATLAB

% Create grouped bar chart
figure;
bar_handle = bar(mean_values, 'grouped');
hold on;

% Apply custom colors to each group
for i = 1:length(bar_handle)
    bar_handle(i).FaceColor = custom_colors(i, :);
end

% Get number of groups and bars per group
[num_groups, num_bars] = size(mean_values);
group_width = min(0.8, num_bars/(num_bars + 1.5)); % Control bar width

% Add error bars
for i = 1:num_bars
    % X positions for each set of bars
    x_positions = bar_handle(i).XEndPoints;
    errorbar(x_positions, mean_values(:, i), std_values(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
    
    % Add significance asterisk if the value is significant
    for j = 1:num_groups
        if significant_values(j, i) == 1
            text(x_positions(j), mean_values(j, i) + std_values(j, i) + 0.02, '*', 'FontSize', 20, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'red');
        end
    end
end

% Set labels and title
xticks(1:num_groups);
xticklabels(categories);
%xlabel('Group');
ylabel('Mean Value');
%title('Mean Values with Standard Deviations and Significance Asterisks (Custom Colors)');
legend(categories, 'Location', 'bestoutside');
ylim([0, 0.38])

hold off;
