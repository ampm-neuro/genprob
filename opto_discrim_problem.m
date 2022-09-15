function [discrims_optoOFF, discrims_optoON, session_files] = opto_discrim_problem(subjects, problem, varargin)
% computes discrimination for multi subjects on first and last day of a
% problem. Makes errorbar barplot


%% input file folder
if nargin == 4
    fp = varargin{1};
else
    fp = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\'];
end



%% first, last
[discrims_optoOFF_first, discrims_optoON_first, session_files_first] = opto_discrim_multisubj(subjects, problem, 1, fp);
[discrims_optoOFF_last, discrims_optoON_last, session_files_last] = opto_discrim_multisubj(subjects, problem, 2, fp);



%% prepare output
discrims_optoOFF = [discrims_optoOFF_first discrims_optoOFF_last];
discrims_optoON = [discrims_optoON_first discrims_optoON_last];
session_files = [session_files_first session_files_last];



%% plot

% remove incomplete problem data
nan_idx = isnan(sum([discrims_optoOFF discrims_optoON],2));
discrims_optoOFF_plot = discrims_optoOFF;
discrims_optoOFF_plot(nan_idx,:) = nan;
discrims_optoON_plot = discrims_optoON;
discrims_optoON_plot(nan_idx,:) = nan;

% edit legend
legend_subjs = cell(1,size(subjects,1));
for isubj = 1:size(subjects,1)
    
    % add subect to legend cell
    legend_subjs(isubj) = subjects(isubj, :);
        
    % annotate subject name if data incomplete
    if nan_idx(isubj)==1
        legend_subjs{isubj} = [legend_subjs{isubj} ' (incomplete)'];
    end
end

% prepare cell for errorbar_barplot input
ebp_input = [{discrims_optoOFF_plot(:,1)} {discrims_optoON_plot(:,1)} {discrims_optoOFF_plot(:,2)} {discrims_optoON_plot(:,2)}];

% plot with legend
figure; hold on
legend_trick(distinguishable_colors(size(subjects,1)), '-')
errorbar_barplot(ebp_input, 1);
legend(legend_subjs, 'location', 'northeastoutside')
title(['Problem ' num2str(problem)])
ylim([-4 6])