function [isesh_hold, isesh_cum, p_dist] = p_dist_rwd_sim(num_sesh, num_round, varargin)
%3 trials at each probability in p_dist
%how logn does it take to converge at p_dist

%plotting input
if nargin==3
    plot_on = varargin{1};
else
    plot_on = 0;
end


%recommended input
%num_sesh = 10;
%num_round = 3;


%x axis
unq_frq = [8000,8151,8305,8462,8621,8784,8950,9119,9291,9466,9645,9827,...
10013,10202,10394,10590,10790,10994,11201,11413,11628,11848,12071,12299,...
12531,12768,13009,13254,13505,13760,14019,14284,14554,14828,15108,15393,...
15684,15980,16282,16589,16902,17221,17546,17877,18215,18559,18909,19266,...
19629,20000];

%{
p_dist = [33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,...
    33,33,33,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,33,33,33,...
    33,33,33,33,33,33,33];
%}

%{
p_dist = [33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,...
    33,34,35,38,44,53,65,77,84,84,77,65,53,44,38,35,34,33,33,33,33,33,...
    33,33,33,33,33,33,33];
    %}

p_dist = [25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,...
    25,25,25,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,25,25,25,...
    25,25,25,25,25,25,25];

num_trl = length(p_dist);

%plot
if plot_on == 1
    h1 = figure;hold on;
end

%sessions
isesh_hold = nan(num_sesh,num_trl);
isesh_cum = nan(num_sesh,num_trl);
for isesh = 1:num_sesh

    %round
    itrl_hold = nan(num_round,num_trl);
    for iround = 1:num_round

        %trials
        for itrl = 1:num_trl

            %current probability
            current_p = p_dist(itrl);
            
            %load if rewarded (pass or fail)
            if current_p >= 100*rand(1)
                if rand(1)<=0.83
                    itrl_hold(iround,itrl) = 1;
                else
                    itrl_hold(iround,itrl) = 0;
                end
            else
                itrl_hold(iround,itrl) = 0;
            end
        end
    
    %smooth individual round
    %itrl_hold(iround,:) = smooth(itrl_hold(iround,:),5);
        
    %plot individual rounds
    if plot_on == 1
        figure(h1)
        subplot(num_sesh*num_round, 3, 1+(iround-1)*3+(isesh-1)*num_round*3)
        plot(unq_frq, itrl_hold(iround,:), 'k-')
        axis([min(unq_frq) max(unq_frq) 0 1]);set(gca,'TickLength',[0, 0]);
        set(gca, 'XScale', 'log')
        box off
        yticks([])
        xticks([])
    end
        
    end
    
    %load session outcomes
    %isesh_hold(isesh,:) = smooth(mean(itrl_hold,1),5);
    %{
    isesh_hold(isesh,:) = mean(itrl_hold,1);
    isesh_cum(isesh,:) = mean(isesh_hold(1:isesh,:),1);
    if isesh > 1       
        
        isesh_gencum_hold = isesh_hold;
        count = 1;
        for i = (isesh-1):-1:1
            count=count+1;
            for i2 = 1:count
                isesh_gencum_hold(i,:) = smooth(isesh_gencum_hold(i,:),5);
            end
        end
        isesh_gencum(isesh,:) = mean([isesh_gencum_hold(1:(isesh-1),:); smooth(isesh_hold(isesh,:),5)']);
    else
        isesh_gencum(isesh,:) = smooth(isesh_hold(isesh,:),5);
    end
    %}
    
    
    isesh_hold(isesh,:) = mean(itrl_hold,1);
    isesh_cum(isesh,:) = mean(isesh_hold(1:isesh,:),1);
    if isesh > 1       
        isesh_gencum(isesh,:) = mean([smooth(isesh_gencum(isesh-1,:),5)'; smooth(isesh_hold(isesh,:),5)']);
    else
        isesh_gencum(isesh,:) = smooth(isesh_hold(isesh,:),5)';
    end

%plot session means and cumulative means
if plot_on == 1
    figure(h1)
    subplot(num_sesh, 3, 2+(isesh-1)*3); hold on
    plot(unq_frq, p_dist./100, '-', 'linewidth', 2.5, 'color', [.8 .8 .8])
    plot(unq_frq, isesh_hold(isesh,:), 'k-')
    axis([min(unq_frq) max(unq_frq) 0 1]);set(gca,'TickLength',[0, 0]);
    set(gca, 'XScale', 'log')
    yticks([0 1])
    xticks([min(unq_frq) max(unq_frq)])
end

if plot_on == 1
    figure(h1)
    subplot(num_sesh, 3, 3+(isesh-1)*3); hold on
    plot(unq_frq, p_dist./100, '-', 'linewidth', 2.5, 'color', [.8 .8 .8])
    plot(unq_frq, isesh_cum(isesh,:), 'k-')
    plot(unq_frq, isesh_gencum(isesh,:), 'r-')
    axis([min(unq_frq) max(unq_frq) 0 1]);set(gca,'TickLength',[0, 0]);
    set(gca, 'XScale', 'log')
    yticks([0 1])
    xticks([min(unq_frq) max(unq_frq)])
end    
    
end