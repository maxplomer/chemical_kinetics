function printgetnu
global nuf nur;
[m,n] = size(nuf);
%%%%%%%%%%%m=53
fid = fopen('odesolver.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETNU(KK,II,NUR,NUF)\n');
fprintf(fid, '      INTEGER KK,II\n');
fprintf(fid, '      INTEGER, dimension(KK,II) :: NUR,NUF\n');
for j=1:n
  for i=1:10:m
    if (i+9)<m
      fprintf(fid, '      NUR(%g:%g,%g) = (/\n',i,i+9,j);
      fprintf(fid, '     *%g,%g,%g,%g,%g,%g,%g,%g,%g,%g/)\n',nur(i,j),nur(i+1,j),nur(i+2,j),nur(i+3,j),nur(i+4,j),nur(i+5,j),nur(i+6,j),nur(i+7,j),nur(i+8,j),nur(i+9,j));
    else
      fprintf(fid, '      NUR(%g:%g,%g) = (/\n',i,m,j);
      fprintf(fid, '     *%g',nur(i,j));
      for k=(i+1):m
        fprintf(fid, ',%g',nur(k,j));
      end
      fprintf(fid, '/)\n');
    end
  end
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end

