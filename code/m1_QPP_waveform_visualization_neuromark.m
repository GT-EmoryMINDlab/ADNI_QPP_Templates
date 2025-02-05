%% ========================================================================
%  QPP Waveform Visualization Script
%  ========================================================================
%  This script loads, processes, and visualizes Quasi-Periodic Patterns 
%  (QPPs) from precomputed .mat files. It focuses on analyzing Neuromark 
%  networks, a set of predefined functional brain networks commonly used 
%  in fMRI studies.
%
%  Key Features:
%  - Loads QPP templates from stored .mat files.
%  - Plots waveform representations of QPPs.
%  - Displays network and subdomain structures using color-coded plots.
%  - Computes correlation matrices for different networks.
%
%  Author: Theodore J. LaGrow
%  Created on: 09/25/2024
%  Last Updated: 02/05/2025
% ========================================================================

%% ========================================================================
%  Clear Workspace and Initialize Environment
% ========================================================================
%  This section clears existing variables, closes any open figures, and 
%  resets the command window for a clean start.

clear; clc; close all

%% ========================================================================
%  Load QPP Templates from .mat Files
% ========================================================================
%  The script scans the directory '../qpp_templates/' for .mat files 
%  containing precomputed QPPs. Each file is loaded into a structured 
%  variable where the filename serves as the field name.

% Define the folder path where QPP templates are stored
folderPath = '../qpp_templates/';

% Get a list of all .mat files in the directory
fileList = dir(fullfile(folderPath, '*.mat'));

% Initialize an empty struct to store QPP data
QPP_templates = struct();

% Loop through each .mat file and load its contents
for i = 1:length(fileList)
    % Extract filename without extension
    [~, fileName, ~] = fileparts(fileList(i).name);
    
    % Load the file
    fileData = load(fullfile(folderPath, fileList(i).name));
    
    % Store the loaded data in a structured variable
    QPP_templates.(fileName) = fileData.(fileName);
end

% Display loaded field names to confirm successful loading
disp(fieldnames(QPP_templates));

%% ========================================================================
%  Select QPP Template for Analysis
% ========================================================================
%  Here, the user selects a specific QPP cohort to analyze.
%  The selected QPP is stored in a cell array for further processing.

cohort = 'ADNI_sNC_QPP'; % Define the QPP cohort of interest

qpp = QPP_templates.(cohort); % Retrieve selected QPP data
QPPs{1,1} = qpp;  % Store it in a cell array for easy indexing

%% ========================================================================
%  Load Neuromark Labels
% ========================================================================
%  The script loads predefined Neuromark network labels, which will be 
%  used to categorize the regions of interest (ROIs) for visualization.

s0_init_neuromark_labels  % Load predefined Neuromark network labels

%% ========================================================================
%  Plot QPP Waveforms for All Networks
% ========================================================================
%  This section generates a time-series plot of QPP waveforms for all 
%  networks, summing over ROIs.

figure;

for icn = 1:size(qpp,1)
    plot(squeeze(sum(QPPs{1,1}(icn,:), 1)) / (size(QPPs{1,1}(icn,:), 1))), hold on,
end

xlabel('Time (s)');
ylabel('QPP value');
ylim([-1.5 1.5]);
title("All Subdomain Networks");

%% ========================================================================
%  Plot QPPs with Network-Specific Coloring
% ========================================================================
%  This section colors each QPP waveform based on its corresponding 
%  functional network using the Networks structure.

figure; 
hold on; % Allow overlaying of multiple plots

% Create placeholder plots for the legend
h = gobjects(numel(Networks), 1);
for i = 1:numel(Networks)
    h(i) = plot(NaN, NaN, 'Color', Networks(i).Color, 'LineWidth', 2); % Placeholder for legend
end

% Loop through each Independent Component Network (ICN) and plot
for i = 1:numel(Networks)
    % Retrieve the ICNs and color for this network
    icn_indices = Networks(i).ICNs;
    color = Networks(i).Color;
    
    % Loop through each ICN and plot with the corresponding network's color
    for icn = icn_indices
        plot(squeeze(sum(QPPs{1,1}(icn,:), 1)) / (size(QPPs{1,1}(icn,:), 1)), ...
             'Color', color, 'LineWidth', 0.5);
    end
end

% Configure the plot
xlabel('Time (s)');
ylabel('QPP value');
ylim([-1.5 1.5]);
title("QPP Waveforms with Network-Specific Coloring");

% Add the legend
legend(h, [Networks.Legend], 'Interpreter', 'latex', 'Location', 'best');

hold off; % Release hold


%% ========================================================================
%  Compute and Display Correlation Matrix for Networks
% ========================================================================
%  This section computes and visualizes the correlation between major 
%  Neuromark networks based on QPP signals.

% Extract QPP waveforms from the Networks structure
networkQPPs = cell2mat(arrayfun(@(x) x.QPP, Networks, 'UniformOutput', false)'); 

% Compute correlation matrix
corrMatrix = corr(networkQPPs');

% Generate the correlation matrix heatmap
figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Label the axes with network names
networkNames = {Networks.Name}; % Extract network names dynamically
numNetworks = numel(Networks);

set(gca, 'XTick', 1:numNetworks, 'XTickLabel', networkNames, ...
         'YTick', 1:numNetworks, 'YTickLabel', networkNames, ...
         'TickLabelInterpreter', 'none');
xtickangle(45);

% Display correlation values in matrix cells
for i = 1:numNetworks
    for j = 1:numNetworks
        text(j, i, sprintf('%.2f', corrMatrix(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end

title('Correlation Matrix of Domain Networks');


%% ========================================================================
%  Compare HC and TN Networks
% ========================================================================
%  This section plots the QPP waveforms for HC and TN networks, highlighting 
%  differences and computing their correlation.

% Retrieve HC and TN network data dynamically
hc_index = find(strcmp({Networks.Name}, 'HC')); % Find HC in structure
tn_index = find(strcmp({Networks.Name}, 'TN')); % Find TN in structure

HC_QPP = Networks(hc_index).QPP; % Get HC QPP waveform
TN_QPP = Networks(tn_index).QPP; % Get TN QPP waveform
HC_Color = Networks(hc_index).Color; % Get HC color
TN_Color = Networks(tn_index).Color; % Get TN color

% Compute correlation between HC and TN waveforms
hc_tn_corr = corr(HC_QPP', TN_QPP');

% Plot HC and TN QPP waveforms
figure;
hold on;
plot(HC_QPP, 'Color', HC_Color, 'LineWidth', 2.0);
plot(TN_QPP, 'Color', TN_Color, 'LineWidth', 2.0);

% Configure plot
xlabel('Time (s)');
ylabel('QPP value');
legend('$HC$', '$TN$', 'Interpreter', 'latex', 'Location', 'best');
ylim([-1 1]);
title(['HC-TN Correlation: ' num2str(hc_tn_corr, '%.2f')]);

hold off;


%% ========================================================================
%  Show Major Subdomains in Neuromark
% ========================================================================
%  This section visualizes QPP waveforms for individual subdomains 
%  within Neuromark networks.

figure;
hold on;

% Loop through each subdomain and plot its QPP waveform
for i = 1:numel(Subdomains)
    plot(Subdomains(i).QPP, 'Color', Subdomains(i).Color, 'LineWidth', 2.0);
end

% Configure plot settings
xlabel('Time (s)');
ylabel('QPP value');
ylim([-1 1]);
title('QPP Waveforms for Neuromark Subdomains');

% Add legend dynamically
legend([Subdomains.Legend], 'Interpreter', 'latex', 'Location', 'best');

hold off;


%% ========================================================================
%  Compute and Display Correlation Matrix for Subdomains
% ========================================================================
%  This section computes and visualizes the correlation between Neuromark 
%  subdomains based on QPP signals.

% Extract QPP waveforms from the Subdomains structure
subdomainQPPs = cell2mat(arrayfun(@(x) x.QPP, Subdomains, 'UniformOutput', false)');

% Compute correlation matrix
corrMatrix = corr(subdomainQPPs');

% Generate the correlation matrix heatmap
figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);

% Label the axes with subdomain names dynamically
subdomainNames = {Subdomains.Name}; % Extract subdomain names
numSubdomains = numel(Subdomains);

set(gca, 'XTick', 1:numSubdomains, 'XTickLabel', subdomainNames, ...
         'YTick', 1:numSubdomains, 'YTickLabel', subdomainNames, ...
         'TickLabelInterpreter', 'none');
xtickangle(45);

% Display correlation values in matrix cells
for i = 1:numSubdomains
    for j = 1:numSubdomains
        text(j, i, sprintf('%.2f', corrMatrix(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end

title('Correlation Matrix of Neuromark Subdomains');


%% ========================================================================
% EOF
% ========================================================================
