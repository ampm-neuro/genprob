function trl_mtx = correct_notone(trl_mtx)
% shift "tone on" event time to after HE for a silent session

trl_mtx(:,7) = trl_mtx(:,10) + (trl_mtx(:,7) - trl_mtx(:,6));
