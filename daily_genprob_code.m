%daily genprob code

%file path
%{
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\pilot';
df = 'pilot';
%}

%
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data';
df = 'repeat_probes\train_consistent';
%}

%{
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data';
df = 'four_tone\consistent_d07_28';
%}

%plot training curves
if strcmp(df, 'pilot')
    [all_stage_learn_mtx, all_stage_learn_cell] = ALL_stage_learn(df, 1:10, 2, 1);
elseif contains(df, 'repeat')
    [all_stage_learn_mtx, all_stage_learn_cell] = ALL_stage_learn(df, 1:7, 2, 1);
else
    [all_stage_learn_mtx, all_stage_learn_cell] = ALL_stage_learn(df, 1:6, 2, 1);
end

%plot fit first probes
%figure; mixed_model_fit(get_file_paths_targeted(fp, {'probe', '_01d', '_01_'}), 3); title('first probe (01d)')
%figure; mixed_model_fit(get_file_paths_targeted(fp, {'probe', '_30d', '_01_'}), 3); title('first probe (30d)')

figure; [fixed_effects_p1, random_effects_p1, beta_pvals_p1, coef_names_p1] ...
    = mixed_model_fit(get_file_paths_targeted([fp '\repeat_probes'], {'postprobe', '_01_'}), 3);
figure; [fixed_effects_p2, random_effects_p2, beta_pvals_p2, coef_names_p2] ...
    = mixed_model_fit(get_file_paths_targeted([fp '\repeat_probes'], {'postprobe', '_02_'}), 3);


%{

%plot first probes
colors = distinguishable_colors(7);
figure; hold on; legend_trick(colors, '-')
for ilearn = 1:7
    plot_all_waits(get_file_paths_targeted_III(fp, {'probe' '_01_'}, {'1d',['learn0' num2str(ilearn)]}, {'30d'}), colors(ilearn,:));
end
legend({'NL01', 'NL02', 'NL03', 'NL04', 'NL05', 'NL06', 'NL07'}, 'location', 'northeastoutside')
title('Probe1 by learning; 01d')
rich_bounds

figure; hold on; legend_trick(colors, '-')
for ilearn = 1:7
    plot_all_waits(get_file_paths_targeted_III(fp, {'probe' '_01_'}, {'30d',['learn0' num2str(ilearn)]}, {'01d'}), colors(ilearn,:));
end
legend({'NL01', 'NL02', 'NL03', 'NL04', 'NL05', 'NL06', 'NL07'}, 'location', 'northeastoutside')
title('Probe1 by learning; 30d')
rich_bounds


%plot second probes
colors = distinguishable_colors(7);
figure; hold on; legend_trick(colors, '-')
for ilearn = 1:7
    plot_all_waits(get_file_paths_targeted_III(fp, {'probe_02'}, {'1d',['learn0' num2str(ilearn)]}, {'30d'}), colors(ilearn,:));
end
legend({'NL01', 'NL02', 'NL03', 'NL04', 'NL05', 'NL06', 'NL07'}, 'location', 'northeastoutside')
title('Probe2 by learning; 01d')
rich_bounds

figure; hold on; legend_trick(colors, '-')
for ilearn = 1:7
    plot_all_waits(get_file_paths_targeted_III(fp, {'probe_02'}, {'30d',['learn0' num2str(ilearn)]}, {'01d'}), colors(ilearn,:));
end
legend({'NL01', 'NL02', 'NL03', 'NL04', 'NL05', 'NL06', 'NL07'}, 'location', 'northeastoutside')
title('Probe2 by learning; 30d')
rich_bounds

%plot third probes
colors = distinguishable_colors(7);
figure; hold on; legend_trick(colors, '-')
for ilearn = 1:7
    plot_all_waits(get_file_paths_targeted_III(fp, {'probe_03'}, {'1d',['learn0' num2str(ilearn)]}, {'30d'}), colors(ilearn,:));
end
legend({'NL01', 'NL02', 'NL03', 'NL04', 'NL05', 'NL06', 'NL07'}, 'location', 'northeastoutside')
title('Probe3 by learning; 01d')
rich_bounds

figure; hold on; legend_trick(colors, '-')
for ilearn = 1:7
    plot_all_waits(get_file_paths_targeted_III(fp, {'probe_03'}, {'30d',['learn0' num2str(ilearn)]}, {'01d'}), colors(ilearn,:));
end
legend({'NL01', 'NL02', 'NL03', 'NL04', 'NL05', 'NL06', 'NL07'}, 'location', 'northeastoutside')
title('Probe3 by learning; 30d')
rich_bounds

%}

%plot new learning curves
openfigs =  findobj('type','figure'); numfigs1 = length(openfigs);
plot_newlearning(df, [1:7], [1], [1:20]); 
openfigs =  findobj('type','figure'); numfigs2 = length(openfigs);
for i = numfigs1+1 : numfigs2
    figure(i); title('01d')
end
openfigs =  findobj('type','figure'); numfigs1 = length(openfigs);
plot_newlearning(df, [1:7], [30], [1:20]); 
openfigs =  findobj('type','figure'); numfigs2 = length(openfigs);
for i = numfigs1+1 : numfigs2
    figure(i); title('30d')
end

%within/between session learning
%wthn_btwn_learn_line
wthn_btwn_learn_line_totals(df, 1); title('01d')
wthn_btwn_learn_line_totals(df, 30); title('30d')

%testpass
[subj_pass, ~, ~, ~, c_sesh] = ALL_testpass(df);
