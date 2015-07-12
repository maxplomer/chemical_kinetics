clc;clear;
%units: erg/(mol)     erg=g*cm^2/s^2=1e-7 J
T=300:50:3000;
chem=ckinit;
[~,k_names]=ckname(chem);
for i=1:length(T)
    h(:,i)=geth(T(i));
end
for i=1:length(T)
    ckh(:,i)=ckhml(T(i),chem);
end
for i=1:9
    figure
    plot(T,h(i,:)*1e-7)
    hold on
    plot(T,ckh(i,:)*1e-7,'.')
    xlabel('T(K)')
    ylabel('h [J/(mol)]')
    title(sprintf('%s',k_names{i}))
end
