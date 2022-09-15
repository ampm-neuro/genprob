function rti = rich_trl_idx(trl_mtx)
% indexes trl_mtx for rich tones unqfrq41(16:26)

% load frequencies
load('unqfrq41', 'unqfrq41')

% floor trl mtx tones
trl_tones = floor(trl_mtx(:,2));

% index
rti = ismember(trl_tones,unqfrq41(16:26));
