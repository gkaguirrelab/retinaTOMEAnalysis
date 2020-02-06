function regressionPlot(xVal, yVal, plotxlabel, plotylabel, plotTitle,plotFlag,debugFlag)

if plotFlag
    if(exist('debugFlag','var'))
        figure(30)
    else
        figure
    end
    plot(xVal,yVal,'xr');
    hold on
    c = polyfit(xVal,yVal,1);
    plot(xVal,polyval(c,xVal),'--b')

    axis square
    xlabel(plotxlabel);
    ylabel(plotylabel);
    title(plotTitle)
end