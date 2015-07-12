clc;clear;
chem=ckinit;
T=1000;
RU=83145100; %erg/(mol*K)
%reaction number
I=1;
g=getg(T);
[nuf,nur]=getnu;
Kp=exp(   -dot(nur(:,I)-nuf(:,I),g)  /  (RU*T)  )
Kc=Kp/(RU*T)^sum(nur(:,I)-nuf(:,I));



nur=cknu(chem)-cknuf(chem);
nuf=-cknu(chem)+nur;

g=ckgml(T,chem);
Kp=exp(   -dot(nur(:,I)-nuf(:,I),g)  /  (RU*T)  )
Kc=Kp/(RU*T)^sum(nur(:,I)-nuf(:,I));




% Kp =
% 
%     5.100545211523303e-018