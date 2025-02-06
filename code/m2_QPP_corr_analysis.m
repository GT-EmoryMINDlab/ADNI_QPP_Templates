%% ========================================================================
%  QPP Correlation and Statistical Analysis Script (Main Script 2)
%  ========================================================================
%  This script aggregates Quasi-Periodic Pattern (QPP) templates across 
%  different cohorts and computes correlation matrices for each template. 
%  It also performs statistical analyses to compare network-level 
%  variations and visualize differences across templates.
%
%  Key Features:
%  - Loads QPP templates dynamically from a predefined structure.
%  - Computes correlation matrices across all 105 networks in each template.
%  - Visualizes the correlation matrix using color-coded plots.
%  - Compares network correlations between selected QPP templates.
%  - Computes variance across networks and identifies statistically 
%    significant differences.
%  - Performs statistical tests (Kruskal-Wallis, Dunn-Sidak, Mann-Whitney U).
%  - Visualizes the correlation differences across intrinsic connectivity 
%    networks (ICNs).
%
%  Author: Theodore J. LaGrow
%  Created on: 09/25/2024
%  Last Updated: 02/06/2025
% ========================================================================

%% ========================================================================
%  Compute and Display Correlation Matrix Across All QPP Templates
%  ========================================================================
%  This script aggregates QPP templates from different cohorts and 
%  computes a correlation matrix across all 105 networks in each template.

% ========================================================================
%  Load and Organize QPP Templates
% ========================================================================
%  The script retrieves QPP data from the `QPP_templates` structure and 
%  arranges them for correlation computation.

% Define cohort template names
templateNames = {'ADNI_sNC_QPP', 'ADNI_sMCI_QPP', 'ADNI_sDAT_QPP', ...
                 'ADNI_uNC_pre_QPP', 'ADNI_uNC_post_QPP', ...
                 'ADNI_pMCI_pre_QPP', 'ADNI_pMCI_post_QPP'};

% Initialize an empty array to store all network data
allNetworks = [];

% Extract QPP templates dynamically from `QPP_templates` structure
for t = 1:numel(templateNames)
    templateName = templateNames{t};
    
    if isfield(QPP_templates, templateName)
        % Retrieve and transpose QPP data for concatenation
        currentQPP = QPP_templates.(templateName)';
        allNetworks = [allNetworks, currentQPP]; % Concatenate across templates
    else
        warning('QPP template "%s" not found in QPP_templates.', templateName);
    end
end

% ========================================================================
%  Compute Correlation Matrix Across All Templates
% ========================================================================
%  This section calculates the correlation matrix for the aggregated
%  QPP data across all templates.

% Compute the correlation matrix
corrMatrix = corr(allNetworks);

% ========================================================================
%  Plot Correlation Matrix Across Templates
% ========================================================================
%  This section visualizes the correlation matrix, adding gridlines
%  to separate different QPP templates.

figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Add lines between templates to visually separate them
hold on;
templateSize = size(QPP_templates.(templateNames{1}), 1); % Assuming all templates have 105 networks

for t = 1:numel(templateNames) - 1
    linePos = t * templateSize;
    xline(linePos, 'k-', 'LineWidth', 5);
    yline(linePos, 'k-', 'LineWidth', 5);
end
axis square;
hold off;

% Add title
title('Correlation Matrix of All 105 Networks Across Templates');

% Save figure
saveas(gcf, 'plots/m2_qpp_template_correlation_matrix.png'); 



%%

%% ========================================================================
%  Compute and Compare QPP Network Correlations
%  ========================================================================

% ========================================================================
%  Load and Compute Correlation Matrices for QPP Templates
% ========================================================================

% Define the two templates for comparison
templateA = 'ADNI_uNC_pre_QPP';
templateB = 'ADNI_pMCI_post_QPP';

% Retrieve QPP data dynamically
if isfield(QPP_templates, templateA)
    A_QPP = QPP_templates.(templateA)(:, 5:12);
else
    error('Template "%s" not found in QPP_templates.', templateA);
end

if isfield(QPP_templates, templateB)
    B_QPP = QPP_templates.(templateB)(:, 5:12);
else
    error('Template "%s" not found in QPP_templates.', templateB);
end

% Compute correlation matrices
A = corr(A_QPP');
B = corr(B_QPP');
A_minus_B = A - B; % Difference matrix

% ========================================================================
%  Plot Correlation Matrices for QPP Templates
% ========================================================================
figure;
subplot(1,3,1);
imagesc(A);
axis square; axis off;
title(['Correlation: ', templateA]);

subplot(1,3,2);
imagesc(B);
axis square; axis off;
title(['Correlation: ', templateB]);

subplot(1,3,3);
imagesc(A_minus_B);
axis square; axis off;
title('Difference (A - B)');

colormap('jet');

% Save figure
saveas(gcf, 'plots/m2_qpp_correlation_comparison.png'); 


% ========================================================================
%  Compute Variance and Statistical Analysis
% ========================================================================
% Compute variance across networks
network_variances = var(A_minus_B, 0, 2); % Row-wise variance

% Expected variance (null hypothesis)
expected_variance = mean(network_variances);

% Chi-square test per network
df = 104; % Degrees of freedom (N-1 where N=105)
chi_sq_values = (network_variances / expected_variance) * df; % Compute chi-square statistic
p_values = 1 - chi2cdf(chi_sq_values, df); % Right-tailed test

% Find significant networks at p < 0.05
significant_networks = find(p_values < 0.05);

% Display results
disp('Significant Networks:');
disp(significant_networks);

% ========================================================================
%  Plot Network Variance Differences
% ========================================================================

% Define domain boundaries for coloring
domain_boundaries = [13, 25, 36, 54, 68, 90, 105];

% Define subdomain boundaries for vertical markers
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
b = bar(network_variances, 'FaceColor', 'flat'); 
hold on;

% Assign colors to bars based on domain structure
for i = 1:numel(Networks)
    idx = Networks(i).ICNs; % Get indices for this domain
    b.CData(idx, :) = repmat(custom_colors(i, :), length(idx), 1); % Assign color
end

% Highlight significant networks
scatter(significant_networks, network_variances(significant_networks), 75, 'r', 'filled'); 

% Add vertical black lines at subdomain boundaries
for i = 1:length(subdomain_boundaries)
    xline(subdomain_boundaries(i) + 0.5, 'k', 'LineWidth', 1.5);
end

% Formatting
axis square;
xlim([0, 105]);
ylim([0, max(network_variances) + 0.1]);

xlabel('Network Index');
ylabel('Variance');
title('Variance Differences Across Networks');

hold off;

% Save figure
saveas(gcf, 'plots/m2_network_variance_comparison.png');



%%


%% ========================================================================
%  Compute and Display Correlation Matrix for QPP Templates
%  ========================================================================

% ========================================================================
%  Load and Organize QPP Templates
% ========================================================================
%  The script retrieves QPP data from `QPP_templates` and ensures availability.

% Define the ordered list of template names
templateNames = {'ADNI_sNC_QPP', 'ADNI_sMCI_QPP', 'ADNI_sDAT_QPP', ...
                 'ADNI_uNC_pre_QPP', 'ADNI_uNC_post_QPP', ...
                 'ADNI_pMCI_pre_QPP', 'ADNI_pMCI_post_QPP'};

numTemplates = numel(templateNames);  % Number of templates

% Initialize cell array to store template data
templates = cell(1, numTemplates);

% Extract templates dynamically
for i = 1:numTemplates
    templateName = templateNames{i};
    if isfield(QPP_templates, templateName)
        templates{i} = QPP_templates.(templateName); % Store the matrix
    else
        error('QPP template "%s" not found in QPP_templates.', templateName);
    end
end

% ========================================================================
%  Compute Correlation Matrix for QPP Templates
% ========================================================================
%  This section calculates the correlation matrix across all QPP templates.

% Initialize the correlation matrix
corrMatrix = zeros(numTemplates, numTemplates);

% Compute correlations
for i = 1:numTemplates
    for j = 1:numTemplates
        % Compute the average absolute correlation between template networks
        corrMatrix(i, j) = mean(mean(abs(corr(templates{i}', templates{j}'))));
    end
end

% Display the correlation matrix
disp('Correlation Matrix:');
disp(corrMatrix);

% ========================================================================
%  Compute Self-Correlation Differences
% ========================================================================
%  This section computes self-correlation minus present value differences.

% Extract self-correlations from the diagonal
selfCorr = diag(corrMatrix);

% Compute the matrix where self-correlation is subtracted from each row
selfCorrMatrix = repmat(selfCorr, 1, numTemplates) - corrMatrix;

% ========================================================================
%  Visualize Correlation Matrices
% ========================================================================
figure;

% First subplot: Regular correlation matrix
subplot(1, 2, 1);
imagesc(corrMatrix);
colormap(jet); 
colorbar;
title('Correlation Matrix');
xlabel('Templates');
ylabel('Templates');
axis square;
set(gca, 'XTick', 1:numTemplates, 'XTickLabel', templateNames, ...
         'YTick', 1:numTemplates, 'YTickLabel', templateNames);
xtickangle(45);
caxis([0, 1]); % Set color scale

% Second subplot: Self-correlation minus present value
subplot(1, 2, 2);
imagesc(selfCorrMatrix);
colormap("autumn");
colorbar;
title('Self-Corr Minus Present Value');
xlabel('Templates');
ylabel('Templates');
axis square;
set(gca, 'XTick', 1:numTemplates, 'XTickLabel', templateNames, ...
         'YTick', 1:numTemplates, 'YTickLabel', templateNames);
xtickangle(45);
caxis([0, 0.3]); % Set color scale

% Save figure
saveas(gcf, 'plots/m2_qpp_template_correlation_analysis.png');




%% ========================================================================
%  Compute and Display Full Correlation Matrices for QPP Templates
%  ========================================================================
%  This script computes full correlation matrices between QPP templates, 
%  calculates self-correlation differences, and visualizes the results.
%

% ========================================================================
%  Load and Organize QPP Templates
% ========================================================================
%  The script retrieves QPP data from `QPP_templates` and ensures availability.

% Define the ordered list of template names
templateNames = {'ADNI_sNC_QPP', 'ADNI_sMCI_QPP', 'ADNI_sDAT_QPP', ...
                 'ADNI_uNC_pre_QPP', 'ADNI_uNC_post_QPP', ...
                 'ADNI_pMCI_pre_QPP', 'ADNI_pMCI_post_QPP'};

numTemplates = numel(templateNames);  % Number of templates

% Initialize cell array to store QPP data
templates = cell(1, numTemplates);

% Extract templates dynamically
for i = 1:numTemplates
    templateName = templateNames{i};
    if isfield(QPP_templates, templateName)
        templates{i} = QPP_templates.(templateName); % Store the matrix
    else
        error('QPP template "%s" not found in QPP_templates.', templateName);
    end
end

% ========================================================================
%  Compute Full and Self Correlation Matrices
% ========================================================================
%  This section calculates self-correlation and full correlation matrices
%  between each pair of QPP templates.

% Initialize storage for correlation matrices
fullCorrMatrices = cell(numTemplates, numTemplates);
selfCorrMatrices = cell(numTemplates, 1);

% Compute correlation matrices
for i = 1:numTemplates
    % Self-correlation for each template
    selfCorrMatrices{i} = corr(templates{i}(:,5:12)', templates{i}(:,5:12)');
    
    for j = 1:numTemplates
        % Full correlation between different templates
        fullCorrMatrices{i, j} = corr(templates{i}(:,5:12)', templates{j}(:,5:12)');
    end
end

% ========================================================================
%  Compute Self-Correlation Difference Matrices
% ========================================================================
%  This section calculates the difference between self-correlation and
%  inter-template correlation matrices.

% Initialize storage for difference matrices
diffMatrices = cell(numTemplates, numTemplates);

% Compute difference from the self-correlation matrix
for i = 1:numTemplates
    for j = 1:numTemplates
        diffMatrices{i, j} = selfCorrMatrices{i} - fullCorrMatrices{i, j};
    end
end

% ========================================================================
%  Visualize Correlation Matrices
% ========================================================================

% Example: Compare sNC vs sMCI (Index 1 vs Index 2)
selected_i = 1; % Index for first template (e.g., sNC)
selected_j = 2; % Index for second template (e.g., sMCI)

figure;

% First subplot: Full correlation matrix between selected templates
subplot(1, 2, 1);
imagesc(fullCorrMatrices{selected_i, selected_j});
colormap(jet);
colorbar;
title(['Full Correlation Matrix (', templateNames{selected_i}, ' vs ', templateNames{selected_j}, ')']);
xlabel('Timepoints/Features');
ylabel('Timepoints/Features');
axis square;

% Second subplot: Difference matrix (Self correlation - full correlation)
subplot(1, 2, 2);
imagesc(diffMatrices{selected_i, selected_j});
colormap(jet);
colorbar;
title(['Difference Matrix (Self Corr - ', templateNames{selected_i}, ' vs ', templateNames{selected_j}, ')']);
xlabel('Timepoints/Features');
ylabel('Timepoints/Features');
axis square;

% Save figure
saveas(gcf, 'plots/m2_qpp_full_correlation_comparison.png');




%% ========================================================================
%  Visualize Difference in Correlation Matrices
% ========================================================================

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

saveas(gcf, 'plots/m2_qpp_difference_correlation_comparison.png'); 


%% ========================================================================
%  Statistical Test
% ========================================================================

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
saveas(gcf, 'plots/m2_boxplot_means_visualization.png'); 


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


%% ========================================================================
%  Visualization of ICN diffences 
% ========================================================================

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
xlabel('Group');
ylabel('Mean Value');
title('Mean Values with Standard Deviations and Significance Asterisks');
legend(categories, 'Location', 'bestoutside');
ylim([0, 0.38])
saveas(gcf, 'plots/m2_qpp_correlation_ICN_significance.png'); 


hold off;

%% ========================================================================
% EOF
% ========================================================================
