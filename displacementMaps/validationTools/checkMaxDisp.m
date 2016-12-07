function checkMaxDisp(radMM,smpPerMM,sectorAngle)

[line,pol,ecc] = createGrid(radMM,smpPerMM);
display('Grid done')
dispMap = makeMap(radMM,smpPerMM,sectorAngle);
display('disp map done')
figure
sz = 25;
scatter(ecc(:),dispMap(:),sz,pol(:),'filled')
axis([0 3 0 3])
hline =refline(1,0)
hline.Color = 'k';
hLine.LineWidth=1.0;
xlabel('Retinal Eccentricity (mm)')
ylabel('Displacement (mm)')
end