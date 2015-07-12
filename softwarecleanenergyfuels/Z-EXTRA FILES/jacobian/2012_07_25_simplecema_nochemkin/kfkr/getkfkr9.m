function [fwdk,revk]=getkfkr9(P,T,X)
%reaction number

I=9;


RU=83145100; %erg/(mol*K)
[A,B,E]=getabe; %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
[nuf,nur]=getnu; 
%third body efficiences
%H2/2.5/ H2O/12/
a=ones(9,1);a(1)=2;a(3)=11;a(2)=0.78;
c=X*P/(RU*T);

%lindemann approach
k0=(6.366E+20)*T^(-1.72)*exp(-(5.248E+02)*41840000/(RU*T));
kinf=A(I)*T^B(I)*exp(-E(I)*41840000/(RU*T));
Pr=(k0*dot(a,c))/kinf;

%Troe form
a=0.8;%alpha
T3=1E-30;%T***
T1=1E+30;%T*
Fcent=(1-a)*exp(-T/T3)+a*exp(-T/T1);   %T** not included
c = -0.4 - 0.67*log(Fcent);
n = 0.75 - 1.27*log(Fcent);
d = 0.14;
logF=log(Fcent)*(1+(  (log(Pr)+c)/(n-d*(log(Pr)+c))  )^2)^-1;
F=exp(logF);
kf=kinf*(Pr/(1+Pr))*F;

%molar concentration c=X*P/(RU*T)
findnuf=find(nuf(:,I));
fwdk=kf;
for i=1:length(findnuf)
   fwdk=fwdk*(X(findnuf(i))*P/(RU*T))^(nuf(findnuf(i),I));
end


%Kp(T)=exp{-[sum i=1:N of (vi''-vi')*mu0,i(T)]/(RU*T)}

g=getg(T);
Kp=exp(   -dot(nur(:,I)-nuf(:,I),g)  /  (RU*T)  );

%kf/kb=Kc=Kp/(RU*T)^[sum i=1:N of vi'' - vi']
PATM=1.01325D6;
Kc=Kp/(RU*T/PATM)^sum(nur(:,I)-nuf(:,I));
kb=kf/Kc;

%molar concentration c=X*P/(RU*T)
findnur=find(nur(:,I));
revk=kb;
for i=1:length(findnur)
   revk=revk*(X(findnur(i))*P/(RU*T))^(nur(findnur(i),I));
end




end