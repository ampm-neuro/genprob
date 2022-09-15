function comp_footprints(orig_folderp, new_folderp)
% plot two sets of footprints, each named the same but in seperate folders

session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];

% original paths 
opaths = get_file_paths_targeted(orig_folderp, {'spatial_footprints'});

% new paths 
npaths = get_file_paths_targeted(new_folderp, {'spatial_footprints'});

% chron sort
opaths = opaths(session_chron_reorder_idx);
npaths = npaths(session_chron_reorder_idx);

% plot each pair
for ipair = 1:size(opaths,1)
    
    % figure
    figure;
    
    % original
    load(opaths{ipair}, 'footprints')
    subplot(1,2,1)
    trace_footprints(footprints)
    
    % new
    load(npaths{ipair}, 'footprints')
    subplot(1,2,2)
    trace_footprints(footprints)
    
    % print
    drawnow
    
end