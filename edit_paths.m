function newpaths = edit_paths(allfp)
% updates filenaming conventions to var standard

if iscell(allfp)
    allfp = allfp(~contains(allfp,'old'));
elseif contains(allfp, 'old')
    newpaths = allfp;
    return
else 
    allfp = {allfp};
end



for ipath = 1:size(allfp,1)

    
    try
    
    if contains(allfp(ipath), '_train')
         cp = allfp{ipath};
         allfp{ipath} = [cp(1:strfind(cp, '_train')) 'lovar' cp(strfind(cp, '_train') + length('_train'): end)];

         %cp = cp
         %ap = allfp{ipath}
         
         movefile( cp, allfp{ipath})
         
    elseif contains(allfp(ipath), '_ctrl')
         cp =  allfp{ipath};
         allfp{ipath} = [cp(1:strfind(cp, '_ctrl')) 'hivar' cp(strfind(cp, '_ctrl') + length('_ctrl'): end)];
         
         %cp = cp
         %ap = allfp{ipath}
         
         movefile( cp, allfp{ipath})
         
    elseif contains(allfp(ipath), '_consistent')
         cp = allfp{ipath};
         allfp{ipath} = [cp(1:strfind(cp, '_consistent')) 'lovar' cp(strfind(cp, '_consistent') + length('_consistent'): end)];
         
         %cp = cp
         %ap = allfp{ipath}
         
         movefile( cp, allfp{ipath})
         
    elseif contains(allfp(ipath), '_scrambled')
         cp =  allfp{ipath};
         allfp{ipath} = [cp(1:strfind(cp, '_scrambled')) 'hivar' cp(strfind(cp, '_scrambled') + length('_scrambled'): end)];
         
         %cp = cp
         %ap = allfp{ipath}
         
         movefile( cp, allfp{ipath})
        
    end
    
    catch
    end
    
    newpaths = allfp{:};
end