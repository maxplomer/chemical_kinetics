function printgetnunet
global nuf nur;
nu=nur-nuf;
[m,n] = size(nuf);
%%%%%%%%%%%m=species
fid = fopen('odesolver.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETNUNET(KK,II,NU)\n');
fprintf(fid, '      INTEGER KK,II\n');
fprintf(fid, '      INTEGER, dimension(KK,II) :: NU\n');
for j=1:n
  for i=1:10:m
    if (i+9)<m
      fprintf(fid, '      NU(%g:%g,%g) = (/\n',i,i+9,j);
      fprintf(fid, '     *%g,%g,%g,%g,%g,%g,%g,%g,%g,%g/)\n',nu(i,j),nu(i+1,j),nu(i+2,j),nu(i+3,j),nu(i+4,j),nu(i+5,j),nu(i+6,j),nu(i+7,j),nu(i+8,j),nu(i+9,j));
    else
      fprintf(fid, '      NU(%g:%g,%g) = (/\n',i,m,j);
      fprintf(fid, '     *%g',nu(i,j));
      for k=(i+1):m
        fprintf(fid, ',%g',nu(k,j));
      end
      fprintf(fid, '/)\n');
    end
  end
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end

