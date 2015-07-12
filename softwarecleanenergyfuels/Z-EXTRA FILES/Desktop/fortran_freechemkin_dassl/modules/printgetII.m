function printgetII
global nur;
[~,II] = size(nur);
fid = fopen('odesolver.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETII (II)\n');
fprintf(fid, '      INTEGER II\n');
fprintf(fid, '      II=%g\n',II);
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');
fclose(fid);
end
