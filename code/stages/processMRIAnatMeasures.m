% Load the biometric and demographic table
ocularMeasuresFileName='~//Dropbox (Aguirre-Brainard Lab)/TOME_subject/TOME-AOSO_SubjectInfo.xlsx';
opts = detectImportOptions(ocularMeasuresFileName);
ocularMeasures = readtable(ocularMeasuresFileName, opts);

idOcularTOME = ocularMeasures.TOME_ID;

age = ocularMeasures.Age;
SE = ocularMeasures.Spherical_Error_average;
axialLength = ocularMeasures.Axial_Length_average;

% Load the OCT results table
octMeasuresFileName ='data/octRGCResultTable.csv';
opts = detectImportOptions(octMeasuresFileName);
octMeasures = readtable(octMeasuresFileName, opts);

% Get the vector of subject IDs
idOctTOME = octMeasures.TOME_ID;

% Load analysis parameters table
anatMeasuresFileName ='data/visualPathwayAnatMeasures.txt';
opts = detectImportOptions(anatMeasuresFileName);
anatMeasures = readtable(anatMeasuresFileName, opts);

% Get the vector of subject IDs
idAnatTOME = anatMeasures.TOME_ID;

% Assemble some anatomical variables
v1Thickness = mean([anatMeasures.lh_lh_v1_noah_template_label_thickness ...
    anatMeasures.rh_rh_v1_noah_template_label_thickness],2);
v1ThicknessScaler = mean([anatMeasures.lh_lh_cortex_label_thickness ...
    anatMeasures.rh_rh_cortex_label_thickness],2);

v1Area = mean([anatMeasures.lh_lh_v1_noah_template_label_area ...
    anatMeasures.rh_rh_v1_noah_template_label_area],2);
areaScaler = mean([anatMeasures.lh_WhiteSurfArea_area ...
    anatMeasures.rh_WhiteSurfArea_area],2);

opticChiasmVolume = anatMeasures.Optic_Chiasm;
periCalcVolume = mean([anatMeasures.wm_lh_pericalcarine ...
    anatMeasures.wm_rh_pericalcarine],2);
volumeScaler = anatMeasures.SupraTentorialVol;

lgnJacobian = anatMeasures.LGNJacobian;

[idAnatInOCT,idOCTInAnat]=matchSubjectID(idAnatTOME, idOctTOME);

corr(opticChiasmVolume,axialLength(idAnatInOcular))

foo = 1;

plot(octMeasures.rgcPCA1(idAnatInOCT(~isnan(idAnatInOCT))),opticChiasmVolume(~isnan(idAnatInOCT)),'xr')

function out = nanForEmpty(in)
    if isempty(in)
        out=nan;
    else
        out=in;
    end
end

function [vecAinB,vecBinA] = matchSubjectID(idVecA, idVecB)

    % Clear out any entries which are not present in both lists
    comboList = [idVecA', idVecB']';
    
    
    vecAinB = cellfun(@(x) nanForEmpty(find(strcmp(x,idVecB))),idVecA);
    vecBinA = cellfun(@(x) nanForEmpty(find(strcmp(x,idVecA))),idVecB);
    
end