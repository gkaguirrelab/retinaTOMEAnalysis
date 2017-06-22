function makeMerdianROIs(params)

% Make 3 ROIs corresponding to the dorsal vertical meridian (dVM), ventral
% vertical meridian (vVM), and horizontal meridan
%
%   Usage:
%       makeMerdianROIs(params)
%
%   Required:
%       params.polFile      = '/path/to/polar/angle/file.nii.gz';
%       params.areasFile    = '/path/to/visual/areas/file.nii.gz';
%       params.hemi         = 'lh' or 'rh';
%       params.dVMFile      = '/path/to/output/dVMfile.nii.gz';
%       params.vVMFile      = '/path/to/output/vVMfile.nii.gz';
%       params.HMFile       = '/path/to/output/HMfile.nii.gz';
%
%   Defaults:
%       params.ROIsize = 15; % degrees polar angle for each ROI
%
%   Written by Andrew S Bock Dec 2016

%% Set defaults
if ~isfield(params,'ROIsize')
    params.ROIsize = 15; % degrees polar angle
end

%% Load in the pRF files
pol     = load_nifti(params.polFile);
areas   = load_nifti(params.areasFile);
%% Get meridian vertices
switch params.hemi
    case 'lh'
        dVM     = find(pol.vol > ((pi/2) - deg2rad(params.ROIsize/2)) & abs(areas.vol) < 3);
        vVM     = find(pol.vol < ((-pi/2) + deg2rad(params.ROIsize/2)) & abs(areas.vol) < 3);
        HM      = find(abs(pol.vol) < deg2rad(params.ROIsize/2) & abs(areas.vol) < 2);
    case 'rh'
        dVM     = find(...
            pol.vol < ((pi/2) + deg2rad(params.ROIsize/2)) & ...
            pol.vol > 0 & ...
            abs(areas.vol) < 3);
        vVM     = find(...
            pol.vol > ((-pi/2) - deg2rad(params.ROIsize/2)) & ...
            pol.vol < 0 & ...
            abs(areas.vol) < 3);
        HM      = find(abs(pol.vol) > (pi - deg2rad(params.ROIsize/2)) & abs(areas.vol) < 2);
end
%% Save the meridian files
% dorsal vertical meridian
out             = pol;
out.vol         = zeros(size(out.vol));
out.vol(dVM)    = 1;
save_nifti(out,params.dVMFile);
% ventral vertical meridian
out             = pol;
out.vol         = zeros(size(out.vol));
out.vol(vVM)    = 1;
save_nifti(out,params.vVMFile);
% horizontal meridian
out             = pol;
out.vol         = zeros(size(out.vol));
out.vol(HM)     = 1;
save_nifti(out,params.HMFile);