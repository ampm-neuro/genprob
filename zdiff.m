function zout = zdiff(v1, v2)
%z difference function returns difference between means of two vectors in
%terms of their pooled variance

if length(v1)>1 && length(v2)>1
    zout = (nanmean(v1(:))-nanmean(v2(:)))/mean([stand_dev(v1(:)) stand_dev(v2(:))]);
    
    
    
elseif length(v1)==1 && length(v2)>1
    zout = (v1-nanmean(v2(:)))/stand_dev(v2(:));
elseif length(v1)>1 && length(v2)==1
    zout = (nanmean(v1(:))-v2)/stand_dev(v1(:));
elseif length(v1)==1 && length(v2)==1
    zout = nan;
else
    error('check input vectors')
end
end