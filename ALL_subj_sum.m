function ALL_subj_sum(training_group)
cur_dir = struct2cell(dir(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\' training_group]));
cur_dir = cur_dir(1,:)';
try % catch missing files
    cur_dir(1:2) = [];
catch
    return
end

%cur_dir
for i = 1:length(cur_dir)
    %cur_dir{i}
    
    try
        subj_summary(cur_dir{i}, training_group)
    catch
        display('no file')
    end
    
end
