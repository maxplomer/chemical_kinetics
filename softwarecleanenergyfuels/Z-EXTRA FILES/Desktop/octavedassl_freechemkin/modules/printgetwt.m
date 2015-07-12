function printgetwt
global species composition;

elementsWT=findelementsWT;

for i=1:length(species)
    WT(i)=dot(composition(i,:),elementsWT);
end

fid = fopen('hardcode/getwt.m', 'wt');
fprintf(fid, 'function [WT]=getwt\n');
fprintf(fid, '%%units: g/mol\n');
fprintf(fid, 'WT=[');
for i=1:length(WT)
    fprintf(fid, '%.6f',WT(i));
    if i<length(WT)
        fprintf(fid, '\n');
    end 
end
fprintf(fid, '];\n');
fprintf(fid, 'end');
fclose(fid);
end