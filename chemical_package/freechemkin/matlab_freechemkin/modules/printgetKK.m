function printgetKK
global species;
KK = length(species);
fid = fopen('hardcode/getKK.m', 'wt');
fprintf(fid, 'function [KK]=getKK\n');
fprintf(fid, '%%number of species\n');
fprintf(fid, 'KK=%5d;\n',KK);
fprintf(fid, 'end');
fclose(fid);
end