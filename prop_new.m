function [props_vect] = prop_new(crm)
% computes proportion of cells in each session that were not active in any
% prior session

% change cell numbers to the number of times seen
for icell = 1:size(crm,1)
    crm(icell,crm(icell,:)>0) = 1:sum(crm(icell,:)>0);
end

% compute proportion of cells that are new
props_vect = nan(1,size(crm,2));
for isesh = 1:size(crm,2)
    props_vect(isesh) = sum(crm(:,isesh)==1)/sum(crm(:,isesh)>0);
end