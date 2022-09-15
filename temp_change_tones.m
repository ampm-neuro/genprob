function [U_new, S_new] = temp_change_tones(U_old, S_old, new_tone_order, sesh_num)

new_tone_order = new_tone_order(sesh_num,:);

old_tone_order = unique(U_old, 'stable')

U_new = nan(size(U_old));
S_new = nan(size(S_old));

for i = 1:length(old_tone_order)
    U_new(U_old==old_tone_order(i)) = new_tone_order(i);
    S_new(U_old==old_tone_order(i)) = mode(S_old(U_old==new_tone_order(i)));
end

csvwrite(['ctl_' num2str(sesh_num) '_S.csv'], S_new);
csvwrite(['ctl_' num2str(sesh_num) '_U.csv'], U_new);
