function h=profilePlot(xVal, yVal, yMean, plotxlabel, plotylabel, plotTitle,plotFlag)

if plotFlag
    h=figure;
    plot(xVal,yVal,'-r');
    hold on
    if ~isempty(yMean)
        plot(xVal,yMean,'-k','LineWidth',4);
    end
    xlabel(plotxlabel);
    ylabel(plotylabel);
    title(plotTitle)
end