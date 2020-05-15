function h=profilePlot(xVal, yVal, yMean, plotxlabel, plotylabel, plotTitle,plotFlag)

if plotFlag
    h=figure;
    if ~isempty(yMean)
        plot(xVal,yVal,'-r');
        hold on
        plot(xVal,yMean,'-k','LineWidth',4);
    else
        plot(xVal,yVal);
    end
    xlabel(plotxlabel);
    ylabel(plotylabel);
    title(plotTitle)
end