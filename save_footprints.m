function save_footprints(folder)
% loads cnmfe_out.mat, computes footprints from variable A, then saves
% footprints as a new mat file.

% load
load([folder '\cnmfe_out.mat'], 'A', 'Cn', 'rev_traces')

% compute footprints
footprints = nan(sum(rev_traces), size(Cn,1), size(Cn,2));

% approved traces
apr_trc = find(rev_traces==1);

% reshape into matrices
for itrace = 1:length(apr_trc)

    % current approved trace
    cat = apr_trc(itrace);
    
    % reshape
    footprints(itrace,:,:) = reshape(A(:,cat), size(Cn,1), size(Cn,2));

    % smooth
    footprints(itrace,:,:) = smooth2a(squeeze(footprints(itrace,:,:)),3);

end

% report
% num_cells = size(footprints,1);

% save
save([folder '\spatial_footprints.mat'], 'footprints');