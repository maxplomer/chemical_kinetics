clc;clear;
%cp=cv + RU
%units: erg/(mol*K)     erg=g*cm^2/s^2=1e-7 J
T=300:50:3000;
chem=ckinit;
[~,k_names]=ckname(chem);
for i=1:length(T)
    cp(:,i)=getcp(T(i));
end
for i=1:length(T)
    ckcp(:,i)=ckcpml(T(i),chem);
end
for i=1:9
    figure
    plot(T,cp(i,:)*1e-7)
    hold on
    plot(T,ckcp(i,:)*1e-7,'.')
    xlabel('T(K)')
    ylabel('cp [J/(mol*K)]')
    title(sprintf('%s',k_names{i}))
end
