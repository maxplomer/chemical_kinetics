function printgetcomposition
global composition;

[m,n] = size(composition);
fid = fopen('mechspecific.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETCOMPOSITION(EE,KK,COMP)\n');
fprintf(fid, '      INTEGER EE, KK\n');
fprintf(fid, '      INTEGER, dimension(EE,KK) :: COMP\n');
for i=1:m
  fprintf(fid, '      COMP(:,%g) = (/%g',i,composition(i,1));
  for j=2:n
    fprintf(fid, ',%g',composition(i,j));
  end
  fprintf(fid, '/)\n');
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end
