%% ========================================================================
%  QPP Projection Script (Main Script 3)
%  ========================================================================
%  This script projects Quasi-Periodic Pattern (QPP) templates onto a 
%  provided fMRI time series dataset and analyzes their occurrence over time.
%
%  Key Features:
%  - Computes correlation between a QPP template and sliding windows of fMRI data.
%  - Tracks occurrences where correlation exceeds a threshold (positive & negative).
%  - Computes occurrences per minute and dwell times between instances.
%  - Generates and saves correlation plots and histograms.
%
%  WARNING: This script requires a time series dataset (D1) to run. If not
%  provided, a random matrix is generated for demonstration purposes.
%
%  Author: Theodore J. LaGrow
%  Created on: 12/06/2024
%  Last Updated: 02/05/2025
% ========================================================================

%% ========================================================================
%  Load or Generate Time Series Data (D1)
% ========================================================================
dataset = 'ADNI_NC_D1';
qpp_proj = 'ADNI_sDAT_QPP';

if ~exist('ADNI_NC_D1', 'var')
    warning('Time series data ADNI_NC_D1 not found. Generating dummy data...');
    ADNI_NC_D1 = randn(105, 1000); % Random matrix (105 ICNs, 1000 time points)
end

if isfield(QPP_templates, qpp_proj)
    qpp = QPP_templates.(qpp_proj)(:, 5:12);  % Extract QPP template
else
    warning('QPP template ADNI_DAT_QPP not found in QPP_templates. Generating dummy template...');
    qpp = randn(105, 8); % Random QPP template (105 ICNs, 8 time points)
end

GrpSn = ADNI_NC_D1;  % Time series data
trrr = 3;  % TR in seconds
corr_thresh = 0.2;  % Correlation threshold

%% ========================================================================
%  Initialize Parameters
% ========================================================================
[~, time_length] = size(GrpSn);
[~, template_length] = size(qpp);

count_above_0_2 = 0;
count_below_neg_0_2 = 0;
correlation_values = zeros(1, time_length - template_length + 1);

%% ========================================================================
%  Compute Sliding Window Correlation
% ========================================================================
for t = 1:(time_length - template_length + 1)
    current_window = GrpSn(:, t:(t + template_length - 1));
    r = corr(current_window(:), qpp(:));
    correlation_values(t) = r; 
    
    if r > corr_thresh
        count_above_0_2 = count_above_0_2 + 1;
    elseif r < -corr_thresh
        count_below_neg_0_2 = count_below_neg_0_2 + 1;
    end
end

%% ========================================================================
%  Compute Occurrences Per Minute
% ========================================================================
total_minutes = (time_length * trrr) / 60;
occurrences_above_0_2_per_min = count_above_0_2 / total_minutes;
occurrences_below_neg_0_2_per_min = count_below_neg_0_2 / total_minutes;

fprintf('Number of occurrences where correlation > 0.2: %d\n', count_above_0_2);
fprintf('Number of occurrences where correlation < -0.2: %d\n', count_below_neg_0_2);
fprintf('Occurrences per minute where correlation > 0.2: %.2f\n', occurrences_above_0_2_per_min);
fprintf('Occurrences per minute where correlation < -0.2: %.2f\n', occurrences_below_neg_0_2_per_min);

%% ========================================================================
%  Plot Correlation Values
% ========================================================================
figure;
plot(correlation_values, 'LineWidth', 1.5);
xlabel('Window Index');
ylabel('Correlation');
title('Correlation Values Across Sliding Windows');
grid on;
saveas(gcf, 'plots/m3_correlation_plot.png');  % Save figure

%% ========================================================================
%  Plot Histogram of Correlation Values
% ========================================================================
figure;
histogram(correlation_values, 'BinWidth', 0.05, 'FaceColor', 'blue', 'EdgeColor', 'black');
xlabel('Correlation Value');
ylabel('Frequency');
title('Histogram of Correlation Values');
grid on;
saveas(gcf, 'plots/m3_correlation_histogram.png');  % Save figure

%% ========================================================================
%  Compute Dwell Times Between Significant Correlation Events
% ========================================================================
above_threshold_indices = find(correlation_values > corr_thresh);
below_threshold_indices = find(correlation_values < -corr_thresh);
all_threshold_indices = sort([above_threshold_indices, below_threshold_indices]);
dwell_times = diff(all_threshold_indices) * trrr;

average_dwell_time = mean(dwell_times, 'omitnan');

fprintf('Average dwell time between instances: %.2f seconds\n', average_dwell_time);
fprintf('Number of dwell time intervals: %d\n', length(dwell_times));

%% ========================================================================
%  Plot Correlation Values With Dwell Time Markers
% ========================================================================
figure;
plot(correlation_values, 'LineWidth', 1.5);
xlabel('Window Index');
ylabel('Correlation');
title('Correlation Values Across Sliding Windows with Dwell Time Markers');
grid on;
saveas(gcf, 'plots/m3_correlation_plot_with_dwell.png');  % Save figure

%% ========================================================================
%  Plot Histogram of Dwell Times
% ========================================================================
figure;
histogram(dwell_times, 'BinWidth', 1, 'FaceColor', 'blue', 'EdgeColor', 'black');
xlabel('Dwell Time Between Instances (seconds)');
ylabel('Frequency');
title('Histogram of Dwell Times Between Instances');
set(gca, 'YScale', 'log'); % Log scale for better visualization
grid on;
saveas(gcf, 'plots/m3_dwell_time_histogram.png');  % Save figure

%% ========================================================================
%  EOF
% ========================================================================