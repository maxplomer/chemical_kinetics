clc;clear;
%units: erg/(mol*K)     erg=g*cm^2/s^2=1e-7 J
T=300:50:3000;
chem=ckinit;
[~,k_names]=ckname(chem);
for i=1:length(T)
    s(:,i)=gets(T(i));
end
for i=1:length(T)
    cks(:,i)=cksml(T(i),chem);
end
for i=1:9
    figure
    plot(T,s(i,:)*1e-7)
    hold on
    plot(T,cks(i,:)*1e-7,'.')
    xlabel('T(K)')
    ylabel('s [J/(mol*K)]')
    title(sprintf('%s',k_names{i}))
end
