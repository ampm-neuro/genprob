function wthn_btwn_learn_line
%plot within and between sesh learning cumulative over every learning stage
%plot sum total learning for each new learning stage

colors = distinguishable_colors(8);
withins_training = nan(10,2);
between_training = nan(10,2);
total_training = nan(10,2);
wb_differences_training = nan(10,2); %mean, sem
for istage = 1:10
    
    istage_hold = num2str(istage);
    if length(istage_hold) < 2
        istage_hold = ['0' istage_hold];
    end
    
    [within_learn, between_learn] = wthn_btwn_learn(['train' istage_hold]); 
    total_learn = nansum([within_learn between_learn],2);

    wb_difference = within_learn - between_learn;
    wb_differences_training(istage,:) = [nanmean(wb_difference) nanstd(wb_difference)./sqrt(sum(~isnan(wb_difference)))];
    
    withins_training(istage,:) = [nanmean(within_learn) nanstd(within_learn)./sqrt(sum(~isnan(within_learn)))];
    between_training(istage,:) = [nanmean(between_learn) nanstd(between_learn)./sqrt(sum(~isnan(between_learn)))];
    total_training(istage,:) = [nanmean(total_learn) nanstd(total_learn)./sqrt(sum(~isnan(total_learn)))];
end

withins_newlearn = nan(10,2);
between_newlearn = nan(10,2);
total_newlearn = nan(10,2);
wb_differences_newlearn = nan(7,2); %mean, sem
for istage = 1:7
    [within_learn, between_learn] = wthn_btwn_learn(['newlearn0' num2str(istage)]);
    total_learn = nansum([within_learn between_learn],2);

    wb_difference = within_learn - between_learn;
    wb_differences_newlearn(istage,:) = [nanmean(wb_difference) nanstd(wb_difference)./sqrt(sum(~isnan(wb_difference)))];
    
    withins_newlearn(istage,:) = [nanmean(within_learn) nanstd(within_learn)./sqrt(sum(~isnan(within_learn)))];
    between_newlearn(istage,:) = [nanmean(between_learn) nanstd(between_learn)./sqrt(sum(~isnan(between_learn)))];
    total_newlearn(istage,:) = [nanmean(total_learn) nanstd(total_learn)./sqrt(sum(~isnan(total_learn)))];
end


%cumulate
%withins_training = nancumsum(withins_training,1);
%between_training = nancumsum(between_training,1);
%total_training = nancumsum(total_training,1); %screws up standard error

%plot withins and betweens
%
figure; hold on
%errorbar(1:10, withins_training(:,1), withins_training(:,2), '-o', 'color', colors(1,:))
%plot(1:10, between_training(:,1), '-', 'color', colors(1,:))
%errorbar(1:10, between_training(:,1), between_training(:,2), '.', 'markersize', 23, 'color', colors(1,:))
plot(1:10, total_training(:,1), '-', 'color', colors(1,:))
errorbar(1:10, total_training(:,1), total_training(:,2), '.', 'markersize', 20, 'color', colors(1,:))
for istage=1:7
    %errorbar(11 + istage*(1/2), withins_newlearn(istage,1), withins_newlearn(istage,2), 'o', 'color', colors(istage+1,:))
    %errorbar(11 + istage*(1/2), between_newlearn(istage,1), between_newlearn(istage,2), '.', 'markersize', 23, 'color', colors(istage+1,:))
    errorbar(11 + istage*(1/2), total_newlearn(istage,1), total_newlearn(istage,2), '.', 'markersize', 20, 'color', colors(istage+1,:))
end
ylabel('Total learning (z)')
%}


xlabel('Training Stage')
set(gca,'TickLength',[0, 0]); box off;
xlim([0 16])
ylim([-5 5])
xticks([1:10 13])
plot(xlim, [1 1].*0, 'k--')
%{
legend({'Training within', '', 'Training between', 'NL01 within', 'NL01 between',...
    'NL02 within', 'NL02 between', 'NL03 within', 'NL03 between', 'NL04 within',...
    'NL04 between', 'NL05 within', 'NL05 between', 'NL06 within', 'NL06 between',...
    'NL07 within', 'NL07 between'}, 'location', 'northeastoutside')
%}
legend({'Training','', 'NL01','NL02','NL03','NL04','NL05','NL06','NL07'}, 'location', 'northeastoutside')



end


function B=nancumsum(A,dim,nmode)
%
if nargin < 3
    nmode = 1;
end
if ~ismember(nmode,1:4)
    error('NANCUMSUM: unacceptable value for nmode parameter.');
end
if nargin < 2 || isempty(dim)
    if ~isscalar(A)
        dim = find(size(A)>1);
        dim = dim(1);
    else
        % For scalar inputs (no nonsingleton dimension)
        dim = 1;
    end
end
% Calculate cumulative sum, depending on selection of nmode
switch nmode
    case 1
        % TREAT NaNs as 0's
        B = A;
        B(B~=B) = 0;
        B = cumsum(B, dim);
    case 2
        % DO NOT INCREMENT, BUT USE NaNs AS PLACEHOLDERS.
        B = nancumsum(A, dim, 1);
        B(A~=A) = NaN;
     case 3
        % RESET sum on NaNs, replacing NaNs with zeros.
        naninds = find(A~=A);
        for ii = 1:numel(naninds)
            B = nancumsum(A, dim, 1);
            A(naninds(ii)) = -B(naninds(ii));
        end
        B = cumsum(A,dim);
    otherwise %case 4
        % RESET sum on NaNs, maintaining NaNs as position holders.
        naninds = find(A~=A);
        B = nancumsum(A,dim,3);
        B(naninds)= NaN;
end

end
