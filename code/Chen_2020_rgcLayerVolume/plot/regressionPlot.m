function h=regressionPlot(xVal, yVal, plotxlabel, plotylabel, plotTitle,plotFlag)

if plotFlag
    if(exist('debugFlag','var'))
        h=figure(30);
    else
        h=figure;
    end
    plot(xVal,yVal,'o','MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on
    c = polyfit(xVal,yVal,1);
    plot(xVal,polyval(c,xVal),'-','Color',[0.5 0.5 0.5])

    axis square
    xlabel(plotxlabel);
    ylabel(plotylabel);
    title(plotTitle)
    box off
end