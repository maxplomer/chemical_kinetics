function printgetwt
global species composition;

elementsWT=findelementsWT;

for i=1:length(species)
    WT(i)=dot(composition(i,:),elementsWT);
end

fid = fopen('mechspecific.f', 'a');

fprintf(fid, '\n      SUBROUTINE GETWT(KK, WT)\n');
fprintf(fid, '      INTEGER KK\n');
fprintf(fid, '      DOUBLE PRECISION, dimension(KK) :: WT\n');
for i=1:length(WT)
  fprintf(fid, '      WT(%g) = %.6f\n',i,WT(i));
end
fprintf(fid, '      RETURN\n');
fprintf(fid, '      END\n');

fclose(fid);
end
