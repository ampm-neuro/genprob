function [doi_m, doi_h] = doi_MevarHivar(doi_h, doi_m)
% overlays output from ALL_probe_v_behavior_fit from both hivar and mevar groups
%{
 [~,~,doi_m] = ALL_probe_vs_behavior_fit('train_mevar_fin');
 [~,~,doi_h] = ALL_probe_vs_behavior_fit('train_hivar_fin');
%}

%
num_subjects = size(doi_h{1},1) + size(doi_m{1},1);
num_problems = size(doi_m{1},2);
subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
stage_num_mtx = repmat(1:6, num_subjects, 1); % stage matrix
model_str = 'var2~predaccuracy+problem+(1|subject)';


%% q1 first day discrim
%allplots(doi_h, doi_m, [1 2], 'q1')


%% q2 first day discrim
%allplots(doi_h, doi_m, [1 3], 'q2')


%% q3 first day discrim
%allplots(doi_h, doi_m, [1 4], 'q3')


%% q4 first day discrim
%allplots(doi_h, doi_m, [1 5], 'q4')
%lme_q1 = lme_function(model_str, [doi_h{1};doi_m{1}], [doi_h{5};doi_m{5}], subj_num_mtx, stage_num_mtx);


%% q5 first day discrim
allplots(doi_h, doi_m, [1 6], 'q5')
%lme_q1 = lme_function(model_str, [doi_h{1};doi_m{1}], [doi_h{5};doi_m{5}], subj_num_mtx, stage_num_mtx);



end


function allplots(h_in, m_in, qidx, title_string)
    figure; hold on;
    %plot2var(h_in{qidx(1)}, h_in{qidx(2)}, winter(size(h_in{qidx(1)},2)))
    %plot2var(m_in{qidx(1)}, m_in{qidx(2)}, autumn(size(h_in{qidx(1)},2)))
    color_buffer = 3;
    hivar_blue = othercolor('Blues9',size(h_in{qidx(1)},2)+color_buffer); hivar_blue = hivar_blue((color_buffer+1):end,:);
    mevar_green = othercolor('Greens9',size(h_in{qidx(1)},2)+color_buffer); mevar_green = mevar_green((color_buffer+1):end,:);
   
    h_in{qidx(1)}(abs(h_in{qidx(1)})>25)=nan;
    h_in{qidx(2)}(abs(h_in{qidx(2)})>25)=nan;
    m_in{qidx(1)}(abs(m_in{qidx(1)})>25)=nan;
    m_in{qidx(2)}(abs(m_in{qidx(2)})>25)=nan;
    
    
    
    plot2var(h_in{qidx(1)}, h_in{qidx(2)}, hivar_blue)
    plot2var(m_in{qidx(1)}, m_in{qidx(2)}, mevar_green)
    
    subplot(1,2,2)
    %[r, p] = fit_line([doi_h{1}(:);doi_m{1}(:)], [doi_h{4}(:);doi_m{4}(:)], 0);
    %title(['r=' num2str(r) ', p=' num2str(p)])
    sgtitle(title_string)
    hold on; [r p ] = fit_line([h_in{qidx(1)}(:); m_in{qidx(1)}(:)],[h_in{qidx(2)}(:); m_in{qidx(2)}(:)], 0)
    xlim([-25 25])
    xlabel('Prediction accuracy (s)')

    figure; hold on
    plot2var_corr([h_in{qidx(1)};m_in{qidx(1)}], [h_in{qidx(2)};m_in{qidx(2)}])
    sgtitle(title_string)
    
    figure; 
    
    
    subplot(1,2,1);hold on
    plot(h_in{qidx(1)}(:,1), h_in{qidx(2)}(:,1), 'ro')
    plot(m_in{qidx(1)}(:,1), m_in{qidx(2)}(:,1), 'bo')
    try
    [r, p] = fit_line([h_in{qidx(1)}(:,1);m_in{qidx(1)}(:,1)], [h_in{qidx(2)}(:,1);m_in{qidx(2)}(:,1)], 0);
     xlim([-25 25])
    plot(xlim, [1 1].*0, 'k--'); plot([1 1].*0, ylim, 'k--')
    title(['first problem; r=' num2str(r) ', p=' num2str(p)])
    catch
    end
    
    
   

    subplot(1,2,2);hold on
    plot(h_in{qidx(1)}(:,end), h_in{qidx(2)}(:,end), 'ro')
    plot(m_in{qidx(1)}(:,end), m_in{qidx(2)}(:,end), 'bo')
    [r, p] = fit_line([h_in{qidx(1)}(:,end);m_in{qidx(1)}(:,end)], [h_in{qidx(2)}(:,end);m_in{qidx(2)}(:,end)], 0);
    xlim([-25 25])
    plot(xlim, [1 1].*0, 'k--'); plot([1 1].*0, ylim, 'k--')
    title(['last problem; r=' num2str(r) ', p=' num2str(p)])
    sgtitle(title_string)
end



%{
% aesthetics
hold on; plot(xlim, [0 0], 'k--')
ylim([-1 1])
tic_vect = [-.8 -.6 -.3 0 .3 .6 .8];
ylim(atanh([-maxr maxr])); 
ylim_hold = ylim;
ylim([ylim_hold(1) ylim_hold(2)+(ylim_hold(2)*.25)])
yticks(atanh(tic_vect)); yticklabels(tic_vect)
%}


% Figure functions
function plot2var(var1, var2, colors)
% vars need to be matrices with (subjects, samples)

var1(abs(var1)>25) = nan;
var2(abs(var2)>25) = nan;

    % colormap
    %colors = bone(10); colors = colors(2:9,:);

    % open figure
    %figure;
    
    % mean centered correlation dot plot
    subplot(1,2,1); hold on
    for iprobe = 1:size(var1,2)

        % mean center
        %
        var1_mc(:, iprobe) = var1(:, iprobe) - nanmean(var1(:, iprobe));
        var2_mc(:, iprobe) = var2(:, iprobe) - nanmean(var2(:, iprobe));
        %}
        plot(var1_mc(:, iprobe), var2_mc(:, iprobe), 'o', 'color', colors(iprobe,:))
    end
    [r, p] = fit_line(var1_mc(:), var2_mc(:), 0);
    set(gca,'TickLength',[0, 0]); box off;
    title(['r=' num2str(r) ', p=' num2str(p)])
    xlabel([remove_underscore(inputname(1)) ', mean centered'])
    ylabel([remove_underscore(inputname(2)) ', mean centered'])
    axis square
    hold on; plot(xlim, [0 0], 'k--')
    hold on; plot([0 0], ylim, 'k--')


    % dot plot with lines
    subplot(1,2,2); hold on
    for isubj = 1:size(var1,1)
        for iprobe = 1:size(var1,2)
            if iprobe < size(var1,2)
                %plot(var1(isubj, [iprobe iprobe+1]), var2(isubj, [iprobe iprobe+1]), '-', 'color', colors(iprobe,:), 'linewidth', 0.5)
            end
            plot(var1(isubj, iprobe), var2(isubj, iprobe), 'o', 'color', colors(iprobe,:))
        end
    end
    for iprobe = 1:size(var1,2)
        if iprobe < size(var1,2)
            %plot(nanmean(var1(:, [iprobe iprobe+1])), nanmean(var2(:, [iprobe iprobe+1])), '-', 'color', colors(iprobe,:), 'linewidth', 6)
        end
        plot(nanmean(var1(:, iprobe)), nanmean(var2(:, iprobe)), '.', 'color', colors(iprobe,:), 'markersize', 60)
    end
    
    set(gca,'TickLength',[0, 0]); box off;
    xlabel(remove_underscore(inputname(1)))
    ylabel(remove_underscore(inputname(2)))
    axis square
    hold on; plot(xlim, [0 0], 'k--')
    hold on; plot([0 0], ylim, 'k--')
    

end

% plot correlations for each probe/problem
function plot2var_corr(var1, var2)

    % if nan in one, nan in both
    nan_idx = isnan(var1) | isnan(var2);
    var1(nan_idx) = nan;
    var2(nan_idx) = nan;

    for iprob = 1:size(var1,2)
        
        if all(isnan(var1(:,iprob))) || all(isnan(var2(:,iprob)))
             continue
        end
        
        subplot(1,size(var1,2),iprob) 
        hold on
        [r,p] = fit_line(var1(:,iprob), var2(:,iprob)); 
        axis square
        plot(xlim, [1 1].*0, 'k--')
        plot([1 1].*0, ylim, 'k--')
        xlabel(remove_underscore(inputname(1)))
        ylabel(remove_underscore(inputname(2)))
        title(['Prob ' num2str(iprob) '; r=' num2str(r) ', p=' num2str(p)])
    end
    set(gcf, 'Position', [680 558 1533 420])
end


function lme = lme_function(model_str, predaccuracy, var2, subj_num_mtx, stage_num_mtx)
    % function testing whether var2 (input 3) predicts var1 (input 2)
    
    
    % if nan in one, nan in both
 %   nan_idx = isnan(var1) | isnan(var2);
 %   var1(nan_idx) = nan_idx;
 %   var2(nan_idx) = nan_idx;
        
    % mixed model
    tbl = table(predaccuracy(:), var2(:), subj_num_mtx(:), stage_num_mtx(:),...
        'VariableNames',{'predaccuracy', 'var2', 'subject', 'problem'});
    tbl.subject = categorical(tbl.subject);
    %tbl.problem = categorical(tbl.problem);
    lme = fitlme(tbl, model_str);
                
end
