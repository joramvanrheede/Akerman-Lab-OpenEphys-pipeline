function pairedlineplot(refvar, compvar, labels, xname, yname)
% function pairedlineplot(refvar, compvar, labels, xname, yname)
% draws a paired line plot between observations made at two different
% conditions of an independent variable
% refvar:   observation set 1 (left); a vector of values
% compvar:  observation set 2 (right); a vector of values
% labels:   labels for refvar and compvar, as two strings in a cell array - {'Ex' 'ample'}
% xname:    x axis label
% yname:    y axis label


plot(ones(length(refvar),1),refvar,'k.','MarkerSize',30)
hold on
plot(ones(length(compvar),1)*2,compvar,'k.','MarkerSize',30)

line([ones(length(refvar),1) ones(length(compvar),1)*2]',[refvar(:) compvar(:)]','Color',[0 0 0],'LineWidth',2)
hold off

xlim([0 3])
set(gca,'XTick',[0 1 2 3])
set(gca,'XTicklabel',{'' labels{:} ''})
set(gca,'FontName','Helvetica','FontSize',16,'FontWeight','bold')
set(gca,'LineWidth',2)
ylabel(yname)
set(gca,'FontName','Helvetica','FontSize',16,'FontWeight','bold')
xlabel(xname)
set(gca,'FontName','Helvetica','FontSize',16,'FontWeight','bold')







