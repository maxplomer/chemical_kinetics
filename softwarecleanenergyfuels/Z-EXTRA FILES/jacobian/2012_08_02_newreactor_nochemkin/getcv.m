function cv=getcv(T)
%cp=cv + RU
%units: erg/(mol*K)     erg=g*cm^2/s^2=1e-7 J
RU=83145100; %erg/(mol*K)
cv=getcp(T)-RU;
cv=cv';


end