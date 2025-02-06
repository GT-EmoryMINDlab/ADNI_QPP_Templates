%% ========================================================================
%  QPP Waveform Visualization Script (Main Script 1)
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
%  Last Updated: 02/06/2025
% ========================================================================

%% ========================================================================
%  Clear Workspace and Initialize Environment
% ========================================================================
%clear; clc; close all

%% ========================================================================
%  Load QPP Templates from .mat Files
% ========================================================================
folderPath = '../qpp_templates/';
fileList = dir(fullfile(folderPath, '*.mat'));

QPP_templates = struct();
for i = 1:length(fileList)
    [~, fileName, ~] = fileparts(fileList(i).name);
    fileData = load(fullfile(folderPath, fileList(i).name));
    QPP_templates.(fileName) = fileData.(fileName);
end
disp(fieldnames(QPP_templates));

%% ========================================================================
%  Select QPP Template for Analysis
% ========================================================================
cohort = 'ADNI_sNC_QPP'; 
qpp = QPP_templates.(cohort); 
QPPs{1,1} = qpp;  

%% ========================================================================
%  Load Neuromark Labels
% ========================================================================
s0_init_neuromark_labels  

%% ========================================================================
%  Plot QPP Waveforms for All Networks
% ========================================================================
figure;
for icn = 1:size(qpp,1)
    plot(squeeze(sum(QPPs{1,1}(icn,:), 1)) / (size(QPPs{1,1}(icn,:), 1))), hold on,
end
xlabel('Time (s)');
ylabel('QPP value');
ylim([-1.5 1.5]);
title("All Subdomain Networks");
saveas(gcf, 'plots/m1_qpp_waveforms.png');  % Save figure

%% ========================================================================
%  Plot QPPs with Network-Specific Coloring
% ========================================================================
figure; hold on;
h = gobjects(numel(Networks), 1);
for i = 1:numel(Networks)
    h(i) = plot(NaN, NaN, 'Color', Networks(i).Color, 'LineWidth', 2);
end
for i = 1:numel(Networks)
    for icn = Networks(i).ICNs
        plot(squeeze(sum(QPPs{1,1}(icn,:), 1)) / (size(QPPs{1,1}(icn,:), 1)), ...
             'Color', Networks(i).Color, 'LineWidth', 0.5);
    end
end
xlabel('Time (s)');
ylabel('QPP value');
ylim([-1.5 1.5]);
title("QPP Waveforms with Network-Specific Coloring");
legend(h, [Networks.Legend], 'Interpreter', 'latex', 'Location', 'best');
saveas(gcf, 'plots/m1_qpp_network_coloring.png');  

%% ========================================================================
%  Compute and Display Correlation Matrix for Networks
% ========================================================================
networkQPPs = cell2mat(arrayfun(@(x) x.QPP, Networks, 'UniformOutput', false)');
corrMatrix = corr(networkQPPs');
figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);
networkNames = {Networks.Name};
numNetworks = numel(Networks);
set(gca, 'XTick', 1:numNetworks, 'XTickLabel', networkNames, ...
         'YTick', 1:numNetworks, 'YTickLabel', networkNames);
xtickangle(45);
for i = 1:numNetworks
    for j = 1:numNetworks
        text(j, i, sprintf('%.2f', corrMatrix(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end
title('Correlation Matrix of Domain Networks');
saveas(gcf, 'plots/m1_network_correlation_matrix.png');  

%% ========================================================================
%  Compare HC and TN Networks
% ========================================================================
hc_index = find(strcmp({Networks.Name}, 'HC')); 
tn_index = find(strcmp({Networks.Name}, 'TN')); 
HC_QPP = Networks(hc_index).QPP; 
TN_QPP = Networks(tn_index).QPP; 
HC_Color = Networks(hc_index).Color; 
TN_Color = Networks(tn_index).Color; 
hc_tn_corr = corr(HC_QPP', TN_QPP');
figure;
hold on;
plot(HC_QPP, 'Color', HC_Color, 'LineWidth', 2.0);
plot(TN_QPP, 'Color', TN_Color, 'LineWidth', 2.0);
xlabel('Time (s)');
ylabel('QPP value');
legend('$HC$', '$TN$', 'Interpreter', 'latex', 'Location', 'best');
ylim([-1 1]);
title(['HC-TN Correlation: ' num2str(hc_tn_corr, '%.2f')]);
saveas(gcf, 'plots/m1_hc_tn_correlation.png');  

%% ========================================================================
%  Show Major Subdomains in Neuromark
% ========================================================================
figure; hold on;
for i = 1:numel(Subdomains)
    plot(Subdomains(i).QPP, 'Color', Subdomains(i).Color, 'LineWidth', 2.0);
end
xlabel('Time (s)');
ylabel('QPP value');
ylim([-1 1]);
title('QPP Waveforms for Neuromark Subdomains');
legend([Subdomains.Legend], 'Interpreter', 'latex', 'Location', 'best');
saveas(gcf, 'plots/m1_subdomain_waveforms.png');  

%% ========================================================================
%  Compute and Display Correlation Matrix for Subdomains
% ========================================================================
subdomainQPPs = cell2mat(arrayfun(@(x) x.QPP, Subdomains, 'UniformOutput', false)');
corrMatrix = corr(subdomainQPPs');
figure;
imagesc(corrMatrix);
colorbar;
colormap('jet');
caxis([-1, 1]);
subdomainNames = {Subdomains.Name};
numSubdomains = numel(Subdomains);
set(gca, 'XTick', 1:numSubdomains, 'XTickLabel', subdomainNames, ...
         'YTick', 1:numSubdomains, 'YTickLabel', subdomainNames);
xtickangle(45);
for i = 1:numSubdomains
    for j = 1:numSubdomains
        text(j, i, sprintf('%.2f', corrMatrix(i,j)), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end
title('Correlation Matrix of Neuromark Subdomains');
saveas(gcf, 'plots/m1_subdomain_correlation_matrix.png');  

%% ========================================================================
% EOF
% ========================================================================