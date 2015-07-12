function printgetthermodata
global thermdata commonT species;

fid = fopen('hardcode/getthermodata.m', 'wt');
fprintf(fid, 'function [thermodata,commonT]=getthermodata\n');

fprintf(fid, 'thermodata=[ ');
for j=1:14
    for i=1:length(species)
        fprintf(fid, '  %s  ',thermdata{i,j});
    end 
    if j<14
        fprintf(fid, '\n     ');
    end
end
fprintf(fid, '];\n');
fprintf(fid, 'commonT=[');
for i=1:length(species)
    fprintf(fid, '%s',commonT{i});
    if i<length(species)
        fprintf(fid, '\n');
    end 
end
fprintf(fid, '];\n');
fprintf(fid, 'end');
fclose(fid);
end