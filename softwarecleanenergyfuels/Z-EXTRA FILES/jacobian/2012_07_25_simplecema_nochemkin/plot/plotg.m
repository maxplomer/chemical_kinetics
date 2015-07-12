clc;clear;
%units: erg/(mol)     erg=g*cm^2/s^2=1e-7 J
T=300:50:3000;
chem=ckinit;
[~,k_names]=ckname(chem);
for i=1:length(T)
    g(:,i)=getg(T(i));
end
for i=1:length(T)
    ckg(:,i)=ckgml(T(i),chem);
end
for i=1:9
    figure
    plot(T,g(i,:)*1e-7)
    hold on
    plot(T,ckg(i,:)*1e-7,'.')
    xlabel('T(K)')
    ylabel('g [J/mol]')
    title(sprintf('%s',k_names{i}))
end
