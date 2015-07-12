function printgetabe
global A B E;

fid = fopen('hardcode/getabe.m', 'wt');
fprintf(fid, 'function [A,B,E]=getabe\n');
%%%%%%%%%
fprintf(fid, 'A=[');
for i=1:length(A)
    fprintf(fid, '%s',A{i});
    if i<length(A)
        fprintf(fid, '\n');
    end 
end
fprintf(fid, '];\n');
%%%%%%%%%
fprintf(fid, 'B=[');
for i=1:length(B)
    fprintf(fid, '%s',B{i});
    if i<length(B)
        fprintf(fid, '\n');
    end 
end
fprintf(fid, '];\n');
%%%%%%%%%
fprintf(fid, 'E=[');
for i=1:length(E)
    fprintf(fid, '%s',E{i});
    if i<length(E)
        fprintf(fid, '\n');
    end 
end
fprintf(fid, '];\n');
%%%%%%%%%
fprintf(fid, 'end');
fclose(fid);
end