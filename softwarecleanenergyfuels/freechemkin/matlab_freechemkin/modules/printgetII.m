function printgetII
global nur;
[~,II] = size(nur);
fid = fopen('hardcode/getII.m', 'wt');
fprintf(fid, 'function [II]=getII\n');
fprintf(fid, '%%number of reactions\n');
fprintf(fid, 'II=%5d;\n',II);
fprintf(fid, 'end');
fclose(fid);
end