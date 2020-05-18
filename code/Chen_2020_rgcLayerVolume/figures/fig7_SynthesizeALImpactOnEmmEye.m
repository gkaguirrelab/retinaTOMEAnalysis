function fig7_SynthesizeALImpactOnEmmEye(scoreExpandedSmoothed,synCoeff,XPos_Degs,saveDir)
%%plot each PC by axial length
h=figure;

%use center profile if we want to see it relative to an emmetropic eye
%(23.58mm) which is in position 16 in "syncoeff"
%centerProfile = scoreExpandedSmoothed(:,1)*synCoeff(16,1)';

for d = 1:6
subplot(3,2,d)

if(d==1)%don't add the center profile if on the first PC
centerProfile = zeros(size(scoreExpandedSmoothed(:,1)*synCoeff(16,1)'));    
else
centerProfile = scoreExpandedSmoothed(:,1)*synCoeff(16,1)';    
end

profileFit = centerProfile + scoreExpandedSmoothed(:,d)*synCoeff(1,d)';
plot(XPos_Degs,profileFit,'-r','LineWidth',1);
hold on
profileFit = centerProfile + scoreExpandedSmoothed(:,d)*synCoeff(16,d)';
plot(XPos_Degs,profileFit,'-g','LineWidth',1);
profileFit = centerProfile + scoreExpandedSmoothed(:,d)*synCoeff(50,d)';
plot(XPos_Degs,profileFit,'-b','LineWidth',1);
xlabel('Eccentricity [deg visual angle]');
xlim([-25 25])
ylim([0 0.007])
setTightFig

end
saveas(h,fullfile(saveDir,'fig7','a.pdf'));

%suptitle('PC components by Axial Lengths: Red=21.79, Green=23.58, Blue=27.57')
