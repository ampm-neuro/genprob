function plot_problemlearning_bysubj
% plots problem discrimination errorbars with each subject contributing one
% d' value



figure; 
for iproblem = [1:6]
    
    subplot(6,1, iproblem); hold on

    % mevar
    plot_all_subj_problem('train_mevar', iproblem, 1, 2)
    plot_all_subj_problem('train_mevar', iproblem, 2, 2)

    % hivar
    plot_all_subj_problem('train_hivar', iproblem, 1, 2)
    plot_all_subj_problem('train_hivar', iproblem, 2, 2)

    title(['Problem ' num2str(iproblem)])
end
drawnow
%sgtitle('hivar')
end
