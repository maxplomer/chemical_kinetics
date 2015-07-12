function printgetspecies
global species;

fid = fopen('hardcode/getspecies.m', 'wt');
fprintf(fid, 'function [species]=getspecies\n');

fprintf(fid, 'species={');
for i=1:length(species)
    fprintf(fid, '''%s''',species{i});
    if i<length(species)
        fprintf(fid, '\n');
    end 
end
fprintf(fid, '};\n');
fprintf(fid, 'end');
fclose(fid);
end