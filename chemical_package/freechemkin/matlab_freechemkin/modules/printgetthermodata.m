function printgetthermodata
global thermdata commonT species;

fid = fopen('hardcode/getthermodata.m', 'wt');
fprintf(fid, 'function [thermodata,commonT]=getthermodata\n');
for i=1:length(species)
    fprintf(fid, '%% %s\n',species{i});
    fprintf(fid, 'thermodata(:,%g)=[ ',i);
    for j=1:5
        fprintf(fid, ' %s ',thermdata{i,j});
    end
    fprintf(fid, ' ... \n');
    for j=6:10
        fprintf(fid, ' %s ',thermdata{i,j});
    end
    fprintf(fid, ' ... \n');
    for j=11:14
        fprintf(fid, ' %s ',thermdata{i,j});
    end
    fprintf(fid, '];\n');
end
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