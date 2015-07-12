function printgetelements
global elements;

fid = fopen('hardcode/getelements.m', 'wt');
fprintf(fid, 'function [elements]=getelements\n');

fprintf(fid, 'elements={');
for i=1:length(elements)
    fprintf(fid, '''%s''',elements{i});
    if i<length(elements)
        fprintf(fid, '\n');
    end 
end
fprintf(fid, '};\n');
fprintf(fid, 'end');
fclose(fid);
end