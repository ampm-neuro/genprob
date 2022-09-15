function sort_file_medass(folder_path)
% sorts medass data files into condition folders and then runs
% file_medass_data on each condition folder

%identify folder contents
immediate_files_and_folders = dir(folder_path);

%delete shit folder things
delete_idx = [];
for iitem = 1:size(immediate_files_and_folders,1)
    if strcmp(immediate_files_and_folders(iitem).name, '.') || ...
            strcmp(immediate_files_and_folders(iitem).name, '..') || ...
            strcmp(immediate_files_and_folders(iitem).name, 'old')
        delete_idx = [delete_idx; iitem];
    end
end
immediate_files_and_folders(delete_idx) = [];

% get path strings
immediate_files_and_folders_hold = struct2cell(immediate_files_and_folders);
immediate_files_and_folders = cell(size(immediate_files_and_folders_hold,2),1);
for icell = 1:size(immediate_files_and_folders_hold,2)
    immediate_files_and_folders{icell} = [immediate_files_and_folders_hold{2,icell}...
        '\' immediate_files_and_folders_hold{1,icell}];
end

% files only
immediate_files = immediate_files_and_folders(isfile(immediate_files_and_folders));
imm_file_names = regexp(immediate_files,filesep,'split');
for im = 1:size(imm_file_names,1)
    imm_file_names{im} = imm_file_names{im}{end};
end

% delete immediate files that are already filed
existing_files = get_file_paths_all(folder_path);
existing_files = existing_files(~contains(existing_files, 'old'));
delete_idx = false(size(immediate_files,1),1);
for iif = 1:size(immediate_files,1)
    if sum(contains(existing_files,imm_file_names{iif}))>1        
        disp('already filed')
        delete(immediate_files{iif}); % delete file
        delete_idx(iif) = true;
    end
end
immediate_files(delete_idx) = [];

% catch empty folder
if isempty(immediate_files)
    disp('no files')
    %return
end

%% group assignments 

% cage sort
cage_sort = cell(6,2);

% train medium variance 2021
cage_sort{1,1} = 'train_mevar_2021';
cage_sort{1,2} = [{'682611'}; {'682612'}; {'682614'}; {'682947'}];

% train high variance 2021
cage_sort{2,1} = 'train_hivar_2021';
cage_sort{2,2} = [{'682613'}; {'682946'}; {'682952'}; {'682953'}; {'689834m3'};...
    {'689834m4'}; {'690330m3'}; {'690330m4'}; {'690331m2'}; {'690331m3'}; ...
    {'691359m3'}; {'691359m4'}; {'691360m2'}; {'691360m3'}; {'691360m4'}; ...
    {'699437m3'}; {'698342m3'}; {'698342m4'};];

% train medium variance imaging HPC
cage_sort{3,1} = 'train_mevar_imaging_hpc';
cage_sort{3,2} = [{'651049m1'}; {'658648m2'}; {'683470m1'}; {'683470m2'};...
    {'683472m1'}; {'679465m1'}; {'679465m2'}; {'690330m1'}; {'690330m2'};...
    {'691359m1'}; {'691359m2'}; {'691360m1'}; {'699437m4'}; {'699438m2'};...
    {'699438m3'}; {'696944m3'}];

% train high variance imaging HPC
cage_sort{4,1} = 'train_hivar_imaging_hpc';
cage_sort{4,2} = [{'683470m3'}; {'683472m2'}; {'683472m3'}; {'687034m4'};...
    {'690331m1'}; {'699437m1'}; {'699437m2'}; {'699438m1'}];
    
% train medium variance opto HIPPOCAMPUS
cage_sort{5,1} = 'train_mevar_optoExp_hpc';
cage_sort{5,2} = [{'654970'}; {'660468'}; {'660469m1'}; {'670571'}; {'675775m2'};...
    {'682945m1'}; {'682945m4'}; {'682615m2'}; {'682615m4'}; {'698346m1'}; ...
    {'698345m1'}; {'698345m2'}; {'698345m3'}; {'698345m4'}; {'698342m1'};...
    {'698342m2'}; {'698346m1'}; {'698346m2'}; {'698346m3'}; {'698346m4'}];

% train high variance opto HIPPOCAMPUS
cage_sort{6,1} = 'train_hivar_optoExp_hpc';
cage_sort{6,2} = [{'682615m3'}; {'682945m2'}; {'682945m3'}];

% cfos medium variance
cage_sort{7,1} = 'train_mevar_cfos';
cage_sort{7,2} = [{'687927m1'}; {'687927m3'}; {'687928m1'}; {'687928m3'}; {'687929m1'};...
    {'687929m3'}; {'687930m1'}; {'687930m3'}; {'689831m1'}; {'689831m3'}; {'689832m1'};...
    {'689832m3'}; {'689833m1'}; {'689833m3'}; {'689834m1'}];

% cfos high variance
cage_sort{8,1} = 'train_hivar_cfos';
cage_sort{8,2} = [{'687927m2'}; {'687927m4'}; {'687928m2'}; {'687928m4'}; {'687929m2'};...
    {'687929m4'}; {'687930m2'}; {'687930m4'}; {'689831m2'}; {'689831m4'}; {'689832m2'};...
    {'689832m4'}; {'689833m2'}; {'689833m4'}; {'689834m2'}];

% mevar six context
cage_sort{9,1} = 'train_mevar_6ctx';
cage_sort{9,2} = [{'704725m1'}; {'704725m2'}; {'704725m3'}; {'704725m4'};...
    {'704726m1'}; {'704726m2'}; {'704726m3'}; {'704726m4'}];



%% graveyard
%{


% train low variance
cage_sort{1,1} = 'train_lovar';
cage_sort{1,2} = [{'666666'}];

% train medium variance
cage_sort{2,1} = 'train_mevar';
cage_sort{2,2} = [{'647761'}; {'647762'}; {'664773'}; {'664774'}; ...
    {'660469m2'}; {'660469m3'}; {'651049m2'}; {'658648m3'}; {'663483'}; {'663484'};...
    {'653900m2'}; {'649961m2'}; {'653900m3'}];

% train medium variance new
cage_sort{3,1} = 'train_mevar_new';
cage_sort{3,2} = [{'677765'}; {'677761'}; {'675776'}; {'675782'}];


% train high variance
cage_sort{5,1} = 'train_hivar';
cage_sort{5,2} = [{'670035'}; {'670041'}; {'670042'}; {'670053'}; {'670054'}];

% train high variance new
cage_sort{6,1} = 'train_hivar_new';
cage_sort{6,2} = [{'675783'}];

% train medium variance imaging HPC
cage_sort{8,1} = 'train_mevar_imaging_hpc';
cage_sort{8,2} = [{'647718'}; {'649961m1'}; {'651049m1'}; {'658648m1'}; {'658648m2'}; {'653900m1'}];


% train medium variance imaging ACC
cage_sort{9,1} = 'train_mevar_imaging_acc';
cage_sort{9,2} = [{'666666'}]; 

% train medium variance opto exp ANTERIOR CINGULATE
cage_sort{11,1} = 'train_mevar_optoExp_acc';
cage_sort{11,2} = [{'670570'}; {'670036'}; {'675775m3'}; {'675775m4'}];

% train medium variance opto ctl
cage_sort{12,1} = 'train_mevar_optoCtl';
cage_sort{12,2} = [{'666666'}];


%}

%% move files into condition folders
for ifolder = 1:size(cage_sort,1)
    
    current_folder = cage_sort{ifolder,1};

    % skip if subject cell empty
    if isempty(cage_sort{ifolder,2})
        continue
    end
    
    % ensure condition folder exists 
    [folder_path '\' current_folder]
    if ~exist([folder_path '\' current_folder], 'dir')
        mkdir([folder_path '\' current_folder])
        %disp(['mkdir ' folder_path '\' current_folder])
    end
    
    % move files
    to_be_moved = immediate_files(contains(immediate_files, cellstr(cage_sort{ifolder,2})));
    for ifile = 1:size(to_be_moved,1)
        try
            movefile(to_be_moved{ifile}, [folder_path '\' current_folder])
        catch
            disp('catch')
            to_be_moved{ifile}
            [folder_path '\' current_folder]
        end
    end
    
    %
    % file
    file_medass_data([folder_path '\' current_folder]); 
    
    % load
    load_medass_all([folder_path '\' current_folder]);
    %}
    
end


% delete empty files (optional)
delete_empty_folders(folder_path)





