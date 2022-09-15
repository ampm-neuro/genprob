function subj_summary(subj_id, training_group)
% plots task progress of a single subject
% outputs general statistics

%% find containing folder

% all subject files
subj_files = get_file_paths_targeted(cd, subj_id);
subj_files = subj_files(~contains(subj_files,'old'));
subj_files = subj_files(contains(subj_files,'two_tone'));
%subj_files = subj_files(contains(subj_files,'meta'));
subj_files = subj_files(contains(subj_files,training_group));

% identify containing folder
cont_fold = cell(size(subj_files,1),1);
for ifile = 1:size(subj_files,1)
    cont_fold{ifile} = subj_files{ifile}(1:strfind(subj_files{ifile}, subj_id)-2);
end
[s,~,j]=unique(cont_fold);
cont_fold = s{mode(j)};

%% find unique learning stages

% all subject files
subj_files = get_file_paths_targeted(cont_fold, subj_id);
subj_files = subj_files(~contains(subj_files,'old'));

% identify learning stages
stage_cell = cell(size(subj_files,1),1);
for ifile = 1:size(subj_files,1)
    gen_str_start = max([strfind(subj_files{ifile}, 'gen0') strfind(subj_files{ifile}, 'gen1')]);
    stage_cell{ifile} = subj_files{ifile}(gen_str_start : gen_str_start + 4);
end
stage_cell = unique(stage_cell);
 
 
%% gather stage files
stage_files = cell(1, length(stage_cell));
for istage = 1:length(stage_cell)
    stage_files{istage} = get_file_paths_targeted(cont_fold, {subj_id, stage_cell{istage}});
    stage_files{istage} = stage_files{istage}(~contains(stage_files{istage},'old'));
end

 
%% plot stage progress
figure; hold on



%set(gcf,'Position', [44 521 2428 829])

set(gcf,'Position', [91         322        1783         597])
%subplot_rows = 2; subplot_cols = sum(~contains(stage_cell, {'14','15'}))+2;
subplot_rows = 2; subplot_cols = sum(~contains(stage_cell, {'14','15'}))+3;

% color probe subplot box
%{
subplot(subplot_rows,subplot_cols, (subplot_rows*subplot_cols-2):(subplot_rows*subplot_cols))
axis_hold = axis;
hold on; rectangle('Position', [axis_hold(1) axis_hold(3) axis_hold(2)-axis_hold(1) axis_hold(4)-axis_hold(3)],...
    'FaceColor', 0.9.*[1 1 1], 'EdgeColor', 0.9.*[1 1 1])
%}

probe_sp_ct = 0;
for istage = 1:length(stage_files)

% reminder sessions
if contains(stage_files{istage}(1), {'remind', 'notone'})
    %{
    subplot(subplot_rows,subplot_cols, subplot_cols)
    mean_probe_waits = nan(size(stage_files{istage},1),2);
    xtl = cell(1,size(stage_files{istage},1));
    for isesh = 1:size(stage_files{istage},1)
        load(stage_files{istage}{isesh}, 'trl_mtx')
        mean_probe_waits(isesh,1) = mean(trl_mtx(trl_mtx(:,3)==0,12)+2);
        mean_probe_waits(isesh,2) = std(trl_mtx(trl_mtx(:,3)==0,12)+2)./sqrt(sum(trl_mtx(:,3)==0));
        xtl{isesh} = ['rem' num2str(isesh)];
    end
    errorbar(mean_probe_waits(:,1), mean_probe_waits(:,2))
    set(gca,'TickLength',[0, 0]); box off;
    xlim([0.5 size(mean_probe_waits,1)+0.5])
    xticklabels(xtl)
    ylim([0 30])
    title('mean reminder waits')
    %}

% probe sessions
elseif contains(stage_files{istage}(1), 'probe_')
    probe_sp_ct = probe_sp_ct+1;
    if probe_sp_ct==1
        subplot(subplot_rows,subplot_cols, (subplot_rows*subplot_cols-2):(subplot_rows*subplot_cols))
    end
    hold on   

    % pre and postprobes
    pres = stage_files{istage}(contains(stage_files{istage}, 'preprobe'));
    posts = stage_files{istage}(contains(stage_files{istage}, 'postprobe'));
    probe_paths = [pres; posts];
    probe_colors = parula(8);
    ALL_coefs = nan(8,4);
    warning ('off','all');
    for iprobe = 1:size(probe_paths,1)
        load(probe_paths{iprobe});
        hold on; [~,~,~,ALL_coefs(iprobe,:)] = wait_times_plot(trl_mtx, 3, probe_colors(iprobe,:));
    end
    warning ('on','all');    
    
    cla
    hold on;
    %ALL_coefs(:,3) = ALL_coefs(:,3).*[1 -1 1 -1 1 -1 1 1]';
    ALL_coefs(:,3) = abs(ALL_coefs(:,3));
    plot(ALL_coefs(:,3))
    plot(ALL_coefs(:,4))
    xlim([.5 8.5]);set(gca, 'XScale', 'linear');xticks auto; xticklabels auto
    plot(xlim, [0 0], 'k--')
    legend({'Log', 'Norm'})
    title('Probe sessions')
    
% learning sessions
else
    % plot zdiffs on first row
    subplot(subplot_rows, subplot_cols, istage) 
    plot_discrim_curves_I(stage_files{istage});
    title(['zdiffs ' stage_cell{istage}])
    ylim([-3 5])
    
    % plot wait time curves on second row
    subplot(subplot_rows, subplot_cols, istage+subplot_cols)
    plot_discrim_curves_II(stage_files{istage});
end
     
end
 
% continuous learning plot
subplot(subplot_rows,subplot_cols, subplot_cols-2:subplot_cols)
continuous_learning(9, subj_id)

