function printFEX
global  species nuf;
[~,II] = size(nuf);
KK=length(species);

fid = fopen('odesolver.f', 'a');

fprintf(fid, '\n      SUBROUTINE  FEX (NEQ, T, Y, YDOT)\n');
fprintf(fid, '      INTEGER  NEQ, KK, II\n');
fprintf(fid, '      DOUBLE PRECISION  T, Temp, Y(%g), YDOT(%g)\n',KK+1,KK+1);
fprintf(fid, '      DOUBLE PRECISION  CV(%g), U(%g), C(%g), WDOT(%g)\n',KK,KK,KK,KK);
fprintf(fid, '      KK=%g\n',KK);
fprintf(fid, '      II=%g\n',II);
fprintf(fid, '      Temp=Y(KK+1)\n');
fprintf(fid, '      C=Y(1:KK)\n');
fprintf(fid, '      CALL GETCV(KK,Temp,CV)\n');
fprintf(fid, '      CALL GETU(KK,Temp,U)\n');
fprintf(fid, '      CALL GETWC(KK,II,Temp,C,WDOT)\n');
fprintf(fid, '      YDOT(1:KK) = WDOT\n');
fprintf(fid, '      YDOT(KK+1) = -DOT_PRODUCT(U,WDOT)/DOT_PRODUCT(C,CV)\n');
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');


fclose(fid);
end
