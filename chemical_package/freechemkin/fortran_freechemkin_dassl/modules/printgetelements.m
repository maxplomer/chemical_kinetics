function printgetelements
global elements;

fid = fopen('mechspecific.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETELEMENTS (EE,ELEMENTS)\n');
fprintf(fid, '      INTEGER EE\n');
fprintf(fid, '      CHARACTER (LEN=3), dimension(EE) :: ELEMENTS\n');
for i=1:length(elements)
  fprintf(fid, '      ELEMENTS(%g) = "%s"\n',i,upper(elements{i}));
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end
