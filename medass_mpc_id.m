function [mpc_id, start_date, subject_id] = medass_mpc_id(mpc_path)
%returns medass mpc file that was running when the medass data file at
%mpc_path was generated


%read output (table)
fid = fopen(mpc_path);
medass_output = textscan(fid, '%s%s%s%s%s%s');
%isolate session ID information
end_id = find(strcmp(medass_output{1,1}, 'MSN:'), 1, 'first');
mpc_id = [];
for mi = 1:size(medass_output,2)
    mpc_id = [mpc_id medass_output{1,mi}(end_id)];
end
mpc_id = mpc_id{find(contains(mpc_id, ':'))+1};


%isolate session start date information
end_id = find(strcmp(medass_output{1,1}, 'Start'), 1, 'first');
start_date = [];
for mi = 1:size(medass_output,2)
    start_date = [start_date medass_output{1,mi}(end_id)];
end
start_date = start_date{find(contains(start_date, ':'))+1};


%isolate subject id information
end_id = find(strcmp(medass_output{1,1}, 'Subject:'), 1, 'first');
subject_id = [];
for mi = 1:size(medass_output,2)
    subject_id = [subject_id medass_output{1,mi}(end_id)];
end
subject_id = subject_id{find(contains(subject_id, ':'))+1};