function summarizeAcuityLogMAR(varargin)
% Reports summary values of acuity measures from the Connectome subjects
%
% Syntax:
%  summarizeAcuityLogMAR
%
% Description:
%   The TOME subjects underwent testing of their foveal visual acuity using
%   the ETDRS distance chart. Acuity was recorded in a decimal (20/20)
%   format, including the number of missed letters. Here I implement a
%   procedure to convert these acuity measures to LogMAR units, and report
%   the average acuity between the two eyes.
%
%   This approach is described here:
%
%          Holladay, Jack T. "Proper method for calculating average visual
%          acuity." Journal of refractive surgery 13.4 (1997): 388-391.
%


%% Parse vargin for options passed here
p = inputParser;

% Optional analysis params
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);

%% Parse and check the parameters
p.parse(varargin{:});

% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);

% The lines on the ETDRS chart
decimalLines = {'10','12.5','16','20','25','32','40','50','63','80','100','125','160','200','250','320','400','500','630','800','1000','1250','1600','2000'};

% The number of letters on each line of the chart
totalLetters = 5;

% The table headers for the two eyes
decimalLabels = {'Acuity__20_x__OD','Acuity__20_x__OS'};

for ii = 1:length(subjectTable.AOSO_ID)
    
    % Loop over the two eyes
    for ee = 1:2
        
        % Obtain the decimal acuity
        decimal = subjectTable.(decimalLabels{ee}){ii};
        
        % Do we have an incomplete line?
        if contains(decimal,{'-','+'})
            [parts,matches] = strsplit(decimal,{'-','+'});
            primaryLine = parts{1};
            nLetters = str2double(parts{2});
            switch matches{1}
                case '+'
                    % Find the next best decimal line
                    nextLine = decimalLines{find(strcmp(decimalLines,primaryLine))-1};
                case '-'
                    % Find the next worse decimal line
                    nextLine = decimalLines{find(strcmp(decimalLines,primaryLine))+1};
            end
            
            % The acuity is the proportional distance between the primry
            % and next line expressed as logMAR
            LogMAR(ii,ee) = -log10(str2double(primaryLine)/20) * ((totalLetters-nLetters)/totalLetters) + ...
                -log10(str2double(nextLine)/20) * ((nLetters)/totalLetters);
            
        else
            LogMAR(ii,ee) = -log10(str2double(decimal)/20);
        end
        
        % Uncomment to report the calculation for every eye
        %{
            str = sprintf(['Subject %d, eye %d, decimal = ' decimal ', logMAR = %2.2f, fractional decimal = %2.2f \n'],ii,ee,LogMAR(ii,ee),20*10^(-LogMAR(ii,ee)));
            fprintf(str);
        %}
        
    end
    
end

% Report the correlation of LogMAR acuity between the two eyes
str = sprintf('The correlation of LogMAR acuity between the two eyes across subjects is R = %2.2f \n',corr(LogMAR(:,1),LogMAR(:,2)));
fprintf(str);

str = sprintf('The median LogMAR acuity (±IQR) across subjects (averaged over eyes) is %2.2f ± %2.3f \n',median(mean(LogMAR,2)),iqr(mean(LogMAR,2)));
fprintf(str);

% Report the correlation of acuity with axial length
str = sprintf('The correlation of LogMAR acuity (averaged over eyes) with axial length is %2.2f \n',corr(mean(LogMAR,2),subjectTable.Axial_Length_average));
fprintf(str);

% Report the worst acuity in the population
str = sprintf('The worst acuity (averaged over eyes) was %2.2f logMAR, or %2.2f decimal \n',min(mean(LogMAR,2)),20*10^(-min(mean(LogMAR,2))));
fprintf(str);


end
