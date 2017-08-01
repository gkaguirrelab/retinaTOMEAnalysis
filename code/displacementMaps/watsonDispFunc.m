function H = watsonDispFunc(meridian, sampleBaseDeg, verbose)
%watsonDispFunc -- Returns the displacment estimates for either the nasal or 
% temporal meridian form the Watson 2014 paper.
%
% Inputs:
%   meridian      = a string indicating the meridian
%                   either 'nasal' or 'temporal'
%   sampleBaseDeg = the sample points to be estimated in degrees
%   verbose       = plot option. set to 'full' to plot.
%
% Ouput:
%   H = the displacement estimates in degrees.
%
%  MAB 2017
  

switch meridian
    case 'nasal'
        % Parameters set by Watson 2014
        x = sampleBaseDeg;
        d = 15.111;
        g = 0.77754;
        u = -0.15933;
        b = 1.7463;
        a = 2.4607 ;
        
        % Watson Equation for Dispalcement 
        TOP = g .* exp(-((x-u)./b).^g) .* ((x-u)./b).^((a.*g) - 1);
        H = d.*( TOP /(b.*gamma(a)));
    case 'temporal'
        % Parameters set by Watson 2014
        x = sampleBaseDeg;
        d = 14.904;
        g = 0.91565;
        u = -0.09386;
        b = 2.4598;
        a = 1.8938;
        
        % Watson Equation for Dispalcement  
        TOP = g .* exp(-((x-u)./b).^g) .* ((x-u)./b).^((a.*g) - 1);
        H = d.*( TOP /(b.*gamma(a)));
    otherwise
        warning('Option for meridian unkown: choices "nasal" or "temporal".')
end
        
% Plot Option 
if strcmp(verbose,'full')
    figure
    plot(x,H)
    xlabel('Eccentricity (deg)')
    ylabel('Displacement (deg)') 
end

end