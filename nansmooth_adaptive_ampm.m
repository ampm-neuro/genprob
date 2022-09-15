function [wdw_mean, wdw_std, wdw_se] = nansmooth_adaptive_ampm(vect_in, min_window_size, min_points_in_window)
% if matrix input, it will smooth each row seperately


wdw_mean = nan(size(vect_in));
wdw_std = nan(size(vect_in));
wdw_se = nan(size(vect_in));

if min(size(vect_in))==1

    for ism = 1:length(vect_in)

        % expand smooth window as needed
        for ibnd = min_window_size:ceil(length(vect_in)/2)

            % current smooth window bounds
            lo_bnd = ism-ibnd; if lo_bnd < 1; lo_bnd=1; end
            hi_bnd = ism+ibnd; if hi_bnd > length(vect_in); hi_bnd=length(vect_in); end

            % current window
            cw = vect_in(lo_bnd:hi_bnd);

            % criteria
            criteria = sum(~isnan(cw))>=min_points_in_window;

            if criteria
                wdw_mean(ism) = nanmean(cw);
                wdw_std(ism) = nanstd(cw);
                wdw_se(ism) = wdw_std(ism)./sqrt(sum(~isnan(cw)));
                break
            else
                continue
            end
        end
    end
    
elseif min(size(vect_in))>1 
    
    for ivect = 1:size(vect_in,1)

        for ism = 1:size(vect_in,2)

            % expand smooth window as needed
            for ibnd = min_window_size:ceil(size(vect_in,2)/2)

                % current smooth window bounds
                lo_bnd = ism-ibnd; if lo_bnd < 1; lo_bnd=1; end
                hi_bnd = ism+ibnd; if hi_bnd > size(vect_in,2); hi_bnd=size(vect_in,2); end

                % current window
                cw = vect_in(ivect, lo_bnd:hi_bnd);

                % criteria
                criteria = sum(~isnan(cw))>=min_points_in_window;

                if criteria
                    wdw_mean(ivect, ism) = nanmean(cw);
                    wdw_std(ivect, ism) = nanstd(cw);
                    wdw_se(ivect, ism) = wdw_std(ivect, ism)./sqrt(sum(~isnan(cw)));
                    break
                else
                    continue
                end
            end
        end 
        
        
    end
    
    
else
    error('input vector size error')
end