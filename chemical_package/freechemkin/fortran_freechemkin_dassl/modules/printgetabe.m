function printgetabe
global A B E;

fid = fopen('mechspecific.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETABE (II, ABE)\n');
fprintf(fid, '      INTEGER II\n');
fprintf(fid, '      DOUBLE PRECISION, dimension(3,II) :: ABE\n');
for i=1:length(A)
  fprintf(fid, '      ABE(:,%g) = (/%e,%s,%s/)\n',i,log(str2num(A{i})),B{i},E{i});
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end
