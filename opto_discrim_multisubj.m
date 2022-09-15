function [discrims_optoOFF, discrims_optoON, session_files] = opto_discrim_multisubj(subjects, problem, first_last, varargin)
% find discrimination index (cohen D) for rich an poor tones under optoON 
% and optoOFF conditions by iterating through multiple subjects

%% input file folder
if nargin == 4
    fp = varargin{1};
else
    fp = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\'];
end



%% preallocate
discrims_optoOFF = nan(size(subjects,1),1);
discrims_optoON = nan(size(discrims_optoOFF));
session_files = cell(size(discrims_optoOFF));



%% iteratively compute and load
for isubj = 1:size(subjects,1)
    current_subj = subjects(isubj,:);
    
    % compute and load
    [discrims_optoOFF(isubj), discrims_optoON(isubj), session_files{isubj}]...
        = opto_discrim_singlesubj(current_subj, problem, first_last);
    
end