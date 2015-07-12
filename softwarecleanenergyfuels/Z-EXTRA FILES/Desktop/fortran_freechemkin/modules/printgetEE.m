function printgetEE
global elements;
EE = length(elements);
fid = fopen('odesolver.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETEE (EE)\n');
fprintf(fid, '      INTEGER EE\n');
fprintf(fid, '      EE=%g\n',EE);
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end
