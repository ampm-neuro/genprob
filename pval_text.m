function pval_text(pval, xpos, ypos)
%plot pval

for ipval = 1:length(pval)
    if pval(ipval) < 0.001
        text(xpos(ipval),ypos(ipval),'p < 0.001','FontSize',10)
    elseif pval(ipval) < 0.01
        text(xpos(ipval),ypos(ipval),'p < 0.01','FontSize',10)
    elseif pval(ipval) < 0.05
        text(xpos(ipval),ypos(ipval),'p < 0.05','FontSize',10)
    else
        text(xpos(ipval),ypos(ipval),['p = ' num2str(round(pval(ipval)*100)/100)], 'FontSize',10)
    end
end

