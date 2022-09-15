function image_similarity_to_prior_probe(subject_id, cell_regist_mtx, tses)
% plots the relationship between neural similarity on this trial and the
% average of the prior probe against the wait times on that trial. Seperate
% plots are made for 'poor' and 'rich' tones.



%% compute values

% iterate through problems
all_poor_tone_rvals = [];
all_rich_tone_rvals = [];
all_poor_tone_waits = [];
all_rich_tone_waits = [];
all_rich_tone_prob_nums = [];
problem_idx_poor = [];
problem_idx_rich = [];
for iproblem = 1:6
    
    % compute r vals and wait times
    try
        current_problem = iproblem+8;
        prior_probe = iproblem;
        next_probe = iproblem+1;
    [poor_tone_rvals, rich_tone_rvals, poor_tone_waits, rich_tone_waits, rich_poor_tone_nums] = ...
        image_corr_popCorr_waits(subject_id, current_problem, [prior_probe next_probe], cell_regist_mtx, tses);
        
        % set means to zero
        %poor_tone_rvals = poor_tone_rvals - mean(poor_tone_rvals);
        %rich_tone_rvals = rich_tone_rvals - mean(rich_tone_rvals);
        %poor_tone_waits = poor_tone_waits - mean(poor_tone_waits);
        %rich_tone_waits = rich_tone_waits - mean(rich_tone_waits);
        
        % load
        all_poor_tone_rvals = [all_poor_tone_rvals; poor_tone_rvals];
        all_rich_tone_rvals = [all_rich_tone_rvals; rich_tone_rvals];
        all_poor_tone_waits = [all_poor_tone_waits; poor_tone_waits];
        all_rich_tone_waits = [all_rich_tone_waits; rich_tone_waits];
        all_rich_tone_prob_nums = [all_rich_tone_prob_nums; rich_poor_tone_nums];
        problem_idx_poor = [problem_idx_poor; repmat(iproblem, size(poor_tone_waits))];
        problem_idx_rich = [problem_idx_rich; repmat(iproblem, size(rich_tone_waits))];
        
    catch
        all_rich_tone_prob_nums = [all_rich_tone_prob_nums; [{1} {1}]];
        disp(['error on iproblem ' num2str(iproblem)])
    end
    
    
    
end



%% plot

% colors
colors = parula(71); colors = colors(26:66,:);

figure

% overall fit lines
subplot(1, 2, 1)
[r, p] = fit_line(all_poor_tone_rvals, all_poor_tone_waits);
title(['Poor tone trials; r=' num2str(round(r*100)/100) ';p=' num2str(round(p*1000)/1000)])
subplot(1, 2, 2)
[r, p] = fit_line(all_rich_tone_rvals, all_rich_tone_waits);
title(['Rich tone trials; r=' num2str(round(r*100)/100) ';p=' num2str(round(p*1000)/1000)])

% iterate through problems
for iproblem = 1:6
    
    % problem idx
    pidx_poor = problem_idx_poor == iproblem;
    pidx_rich = problem_idx_rich == iproblem;
    
    % poor
    subplot(1, 2, 1)
    hold on
    current_color = colors(all_rich_tone_prob_nums{iproblem, 2},:);
    plot(all_poor_tone_rvals(pidx_poor), all_poor_tone_waits(pidx_poor), 'o', 'color', current_color(1,:))
    
    % rich
    subplot(1, 2, 2)
    hold on
    current_color = colors(all_rich_tone_prob_nums{iproblem, 1},:);
    plot(all_poor_tone_rvals(pidx_rich), all_poor_tone_waits(pidx_rich), 'o', 'color', current_color(1,:))
    
end