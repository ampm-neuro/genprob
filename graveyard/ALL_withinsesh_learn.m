function ALL_withinsesh_learn(stage, sesh, subj)
%generic function that runs through sessions

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen_richards';


%session selection
%stage = 'gt05'
%sesh = '01'
%subj = 'm4'%{'m3', 'm1'}


%get trial matrices
%
%preallocate
%
trl_mtxs = cell(1);
medcells = cell(1);
mwd = cell(1);
count = 0;
all_diff = [];
all_z = [];

%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %exclude undesired subjects
    if ~contains(current_subj, subj)
        continue
    end
    
    %session files
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name;

        
        %exclude undesired sessions
       
        if ~contains(current_sesh, stage) || ~contains(current_sesh, sesh)
            continue
        else
            
        end
        
        %load data file
        %[folderpath '\' current_subj '\' current_sesh]
        load([folderpath '\' current_subj '\' current_sesh])
        

        %load output
        count = count+1;
        trl_mtxs{count} = trl_mtx;
        medcells{count} = medass_cell;

        
    end
end

%figure
f1 = figure; hold on
f2 = figure; hold on
colors = distinguishable_colors(count);

%number of bins
num_bins = 2;


%for each sesh
for isesh = 1:length(trl_mtxs)  
    
    %window size (number of trials)
    winsz = ceil(size(trl_mtxs{isesh},1)/(num_bins+1));
    
    %preallocate
    rich_all = cell(1);
    poor_all = cell(1);
    
    %for each window
    window_count = 0;
    for iw = 1:winsz:size(trl_mtxs{isesh},1)
        window_count = window_count+1;
        
        %trl mtx for this sesh and window
        if size(trl_mtxs{isesh},1)>=iw+2*winsz
            first_bin = iw;
            last_bin = iw+winsz-1;
        else
            first_bin = iw;
            last_bin = size(trl_mtxs{isesh},1);
        end
        local_trl_mtx = trl_mtxs{isesh}(first_bin:last_bin,:);
        
        %compute wait durations
        [mean_wait_durations, wait_durations, ~, unq_frq, p_dist] = wait_times(local_trl_mtx,medcells{isesh},0);
        
        mwd{window_count} = mean_wait_durations;
        
        min_samps = 1;
        
        rich_all{window_count} = cell2mat(wait_durations(p_dist>0.50));       
            if length(rich_all{window_count}) < min_samps
                rich_all{window_count} = nan; 
                mwd{window_count}(2) = nan;
            end
            
        poor_all{window_count} = cell2mat(wait_durations(p_dist<0.50));
            if length(poor_all{window_count}) < min_samps
                poor_all{window_count} = nan; 
                mwd{window_count}(1) = nan;
            end
               
        if last_bin == size(trl_mtxs{isesh},1)
            break
        end
        
    end
     
    %zdiffs for all windows
    means_z=nan(2,length(rich_all));stds_z=nan(size(means_z));wait_times_out=nan(2,length(rich_all));
    for iz = 1:length(rich_all)
        means_z(:,iz) = [mean(rich_all{iz}); mean(poor_all{iz})];
        stds_z(:,iz) = [std(rich_all{iz}); std(poor_all{iz})];
        zdists(iz) = (means_z(1,iz)-means_z(2,iz))/mean(stds_z(:,iz));
        wait_times_out(:,iz) = [mean(poor_all{iz}) mean(rich_all{iz})];
    end
    
    %plot
    figure(f1)
    plot(1:size(wait_times_out,2), wait_times_out(1,:), '-o', 'color', colors(isesh,:), 'linewidth', 1)
    plot(1:size(wait_times_out,2), wait_times_out(2,:), '--o', 'color', colors(isesh,:), 'linewidth', 1)
    title seconds
    
    %capture data from each subj
    all_diff = [all_diff; wait_times_out(2,:)-wait_times_out(1,:)]; %rich-poor
    all_z = [all_z; zdists]; %rich-poor
    
end

figure(f1)
ylim_hold = ylim;
ylim([0 ylim_hold(2)])
set(gca,'TickLength',[0, 0]);
xlim([.5 length(zdists)+.5])
ylabel('Wait time (s)')
legend({'Poor tone','Rich tone'}, 'location', 'Northeastoutside')

figure(f2)
for i = 1:size(all_z,1)
    plot(1:length(zdists), all_z(i,:), '-o', 'color', 0.80.*[1 1 1]);
end
errorbar(nanmean(all_z,1), nanstd(all_z,[],1)./sqrt(sum(~isnan(all_z),1)), 'k', 'linewidth', 3)
set(gca,'TickLength',[0, 0]);
ylim([-5 5])
xlim([.5 size(all_z,2)+.5])
plot(xlim, [0 0], 'k--')
ylabel('Preference for rich tone (z)')




%plot all diff
figure; hold on
for i = 1:size(all_diff,1)
   plot(all_diff(i,:), '-o', 'color', 0.80.*[1 1 1]) 
end
errorbar(nanmean(all_diff,1), nanstd(all_diff,[],1)./sqrt(sum(~isnan(all_diff),1)), 'k', 'linewidth', 3)
xlim([.5 size(all_diff,2)+.5])
plot(xlim, [0 0], 'k--')
set(gca,'TickLength',[0, 0]);
ylabel('Preference for rich tone (s)')

nanmean(all_diff)
[~,p_z,~, stats_z] = ttest(all_diff(:,1), all_diff(:,end))

