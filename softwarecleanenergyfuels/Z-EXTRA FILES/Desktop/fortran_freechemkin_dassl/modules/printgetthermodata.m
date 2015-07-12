function printgetthermodata
global thermdata commonT species;
td=thermdata;
fid = fopen('odesolver.f', 'a');


fprintf(fid, '\n      SUBROUTINE GETTHERMODATA(KK,N)\n');
fprintf(fid, '      INTEGER KK\n');
fprintf(fid, '      DOUBLE PRECISION, dimension(14,KK) :: N\n');
for i=1:length(species)
  fprintf(fid, '      N(:,%g) = (/\n',i);
  fprintf(fid, '     *%s, %s, %s,\n',td{i,1},td{i,2},td{i,3});
  fprintf(fid, '     *%s, %s, %s,\n',td{i,4},td{i,5},td{i,6});
  fprintf(fid, '     *%s, %s, %s,\n',td{i,7},td{i,8},td{i,9});
  fprintf(fid, '     *%s, %s, %s,\n',td{i,10},td{i,11},td{i,12});
  fprintf(fid, '     *%s, %s/)\n',td{i,13},td{i,14});
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');

fprintf(fid, '\n      SUBROUTINE GETCOMMONT(KK, COMMONT)\n');
fprintf(fid, '      INTEGER KK\n');
fprintf(fid, '      DOUBLE PRECISION, dimension(KK) :: COMMONT\n');
for i=1:length(species)
  fprintf(fid, '      COMMONT(%g) = %s\n',i,commonT{i});
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end
