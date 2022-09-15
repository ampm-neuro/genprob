function sig_asterisks(pval, xpos, ypos)
%plot asterisks at xpos,ypos according to the pval
% * for < 0.05
% ** for < 0.01
% *** for < 0.001


for ipval = 1:length(pval)
    if pval(ipval) < 0.001
        text(xpos(ipval),ypos(ipval),'***','FontSize',14)
    elseif pval(ipval) < 0.01
        text(xpos(ipval),ypos(ipval),'**','FontSize',14)
    elseif pval(ipval) < 0.05
        text(xpos(ipval),ypos(ipval),'*','FontSize',14)
    else
        pval_text(pval(ipval), xpos(ipval),ypos(ipval))
    end
end

