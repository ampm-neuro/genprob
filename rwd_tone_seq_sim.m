function [bar_vect_non_rwd_trls, all_bar_vect] = rwd_tone_seq_sim(rwd1_nrwd_seqence, unqfrq)
%plots bar of subj experience with tones over 11 trials

%figure vects
bar_vect_rwd = zeros(1,length(unqfrq));
bar_vect_nrwd = zeros(1,length(unqfrq));

%form heatmap for each trial
for itrl = 1:size(rwd1_nrwd_seqence,1)
    
    %rewarded tone position
    rt_pos = find(unqfrq==rwd1_nrwd_seqence(itrl,1));
    
    %nonrewarded tone positions
    nrt_pos = ismember(unqfrq, rwd1_nrwd_seqence(itrl,2:end));
    
    %add rewarded tone
    bar_vect_rwd(rt_pos) = bar_vect_rwd(rt_pos) + 1;
    
    %add non-rewarded tone
    bar_vect_nrwd(nrt_pos) = bar_vect_nrwd(nrt_pos) + 1;
    
    %output
    bar_vect_non_rwd_trls(itrl,rt_pos) = 2;
    bar_vect_non_rwd_trls(itrl,nrt_pos) = 1;
    %bar_vect_non_rwd_trls(itrl,rt_pos-4 : rt_pos+4) = bar_vect_non_rwd_trls(itrl,rt_pos-4 : rt_pos+4)+2;
    
end

all_bar_vect = [bar_vect_nrwd, bar_vect_rwd];

%figure
figure; hold on
bar([bar_vect_nrwd; bar_vect_rwd]', 'stacked')
%bar(bar_vect_nrwd)
set(gca,'TickLength',[0, 0]); box off;
xticks(1:2:length(unqfrq))
xlim([.5 length(unqfrq)+.5])
ylim([0 4.5])

figure;
%imagesc(bar_vect_non_rwd_trls)
data_hold = bar_vect_non_rwd_trls; data_hold(data_hold==0) = nan;
[nr,nc] = size(data_hold);
pcolor([data_hold nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(gca,'TickLength',[0, 0]); box off;
xticks(1:length(unqfrq))
xticks(1.5:2:length(unqfrq)+.5)
xticklabels(1:2:length(unqfrq))
yticks(1.5:size(bar_vect_non_rwd_trls,1)+.5)
yticklabels(1:size(bar_vect_non_rwd_trls,1))
xlim([.5 length(unqfrq)+.5])

%{
figure;
%imagesc(bar_vect_non_rwd_trls)
data_hold = probe_hm; data_hold(data_hold==0) = nan;
[nr,nc] = size(data_hold);
pcolor([data_hold nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(gca,'TickLength',[0, 0]); box off;
xticks(1.5:2:size(probe_hm,2)+.5)
xticklabels(1:2:size(probe_hm,2))
yticks(1.5:size(probe_hm,1)+.5)
yticklabels(1:size(probe_hm,1))
%xlim([.5 length(unqfrq)+.5])
%}


%pdf
%{
rt_pos_idx = find(max(tone_stimuli_heatmap)==2);
figure; plot(linspace(1,41,41), normpdf(linspace(1,41,41), mean(rt_pos_idx), std(rt_pos_idx)));
%}