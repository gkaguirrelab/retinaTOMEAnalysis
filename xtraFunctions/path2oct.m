function outPath = path2oct
path = pwd;
pos = strfind(path,'octAnalysisForTOME');
outPath = path(1:pos+17);
end