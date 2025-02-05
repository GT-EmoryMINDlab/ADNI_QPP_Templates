%% ========================================================================
%  Neuromark Network and Subdomain Label Initialization
%  ========================================================================
%  This script defines the domain and subdomain labels for Neuromark-based
%  Quasi-Periodic Pattern (QPP) analysis. It initializes the structured 
%  representations of networks and their independent component networks (ICNs).
%
%  Key Features:
%  - Defines Neuromark networks and their corresponding ICNs.
%  - Assigns a unique color to each network for visualization.
%  - Computes and stores QPP waveforms for each network.
%  - Defines subdomains within each network and assigns ICNs.
%  - Prepares structured data for streamlined analysis and plotting.
%
%  Author: Theodore J. LaGrow
%  Created on: 02/05/2025
%  Last Updated: 02/05/2025
% ========================================================================

%% ========================================================================
%  Define Domain Structure for QPP Analysis
% ========================================================================
%  This section defines the main functional networks in Neuromark and 
%  organizes their respective ICNs, colors, and QPP waveform calculations.

Networks = struct( ...
    'Name',    {'CB', 'VI', 'PL', 'SC', 'SM', 'HC', 'TN'}, ...   % Network Names
    'ICNs',    {1:13, 14:25, 26:36, 37:54, 55:68, 69:90, 91:105}, ... % ICN Indices
    'Color',   {[0, 1, 1], [0, 1, 0], [0.5, 0.5, 0.5], [1, 0, 1], [1, 0, 0], [0, 0, 1], [1, 0.5, 0]}, ... % Colors
    'Legend',  {"$CB$", "$VI$", "$PL$", "$SC$", "$SM$", "$HC$", "$TN$"} ... % Legend Text
);

% Compute QPP waveforms for each network
for i = 1:numel(Networks)
    icn_indices = Networks(i).ICNs; % Extract ICN indices
    Networks(i).QPP = squeeze(sum(QPPs{1,1}(icn_indices,:), 1)) / numel(icn_indices); % Compute QPP waveform
end

%% ========================================================================
%  Define Subdomain Structure for QPP Analysis
% ========================================================================
%  This section organizes subdomains within each Neuromark network,
%  assigns their respective ICNs, and computes their QPP waveforms.

Subdomains = struct( ...
    'Name',    {'CB-CB', 'VI-OT', 'VI-OC', 'PL-PL', 'SC-EH', 'SC-ET', 'SC-BG', 'SM-SM', ...
                'HC-IT', 'HC-TP', 'HC-FR', 'TN-CE', 'TN-DN', 'TN-SN'}, ...   % Subdomain Names
    'ICNs',    {1:13, 14:19, 20:25, 26:36, 37:39, 40:45, 46:54, 55:68, ...
                69:75, 76:80, 81:90, 91:93, 94:101, 102:105}, ...  % ICN Indices
    'Color',   {[0, 1, 1], [1, 0.5, 0], [0.5, 0, 0.5], [0.5, 0.5, 0.5], [1, 0.5, 0.5], ...
                [1, 1, 0], [0, 1, 0], [1, 0, 1], [0, 0.5, 1], [1, 0.75, 0.8], ...
                [0.75, 1, 0.25], [1, 0, 0], [0, 0, 1], [0.25, 0.25, 0.25]}, ... % Colors
    'Legend',  {"$CB-CB$", "$VI-OT$", "$VI-OC$", "$PL-PL$", "$SC-EH$", "$SC-ET$", "$SC-BG$", "$SM-SM$", ...
                "$HC-IT$", "$HC-TP$", "$HC-FR$", "$TN-CE$", "$TN-DN$", "$TN-SN$"} ... % Legend Text
);

% Compute QPP waveforms for each subdomain
for i = 1:numel(Subdomains)
    icn_indices = Subdomains(i).ICNs; % Extract ICN indices
    Subdomains(i).QPP = squeeze(sum(QPPs{1,1}(icn_indices,:), 1)) / numel(icn_indices); % Compute QPP waveform
end

%% ========================================================================
%  EOF
% ========================================================================
