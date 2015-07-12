function printgetspecies
global species;

fid = fopen('odesolver.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETSPECIES (KK,SPECIES)\n');
fprintf(fid, '      INTEGER KK\n');
fprintf(fid, '      CHARACTER (LEN=16), dimension(KK) :: SPECIES\n');
for i=1:length(species)
  fprintf(fid, '      SPECIES(%g) = "%s"\n',i,upper(species{i}));
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end

