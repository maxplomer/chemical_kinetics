function printgetKK
global species;
KK = length(species);
fid = fopen('odesolver.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETKK (KK)\n');
fprintf(fid, '      INTEGER KK\n');
fprintf(fid, '      KK=%g\n',KK);
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end
