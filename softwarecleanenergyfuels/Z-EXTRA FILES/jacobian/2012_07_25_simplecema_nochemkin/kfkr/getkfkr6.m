function [fwdk,revk]=getkfkr6(P,T,X)
%reaction number
I=6;


RU=83145100; %erg/(mol*K)
[A,B,E]=getabe; %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
[nuf,nur]=getnu;                            
kf=A(I)*T^B(I)*exp(-E(I)*41840000/(RU*T));

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


%third body efficiences
%H2/2.5/ H2O/12/
a=ones(9,1);a(1)=2.5;a(3)=12;
c=X*P/(RU*T);

fwdk=fwdk*dot(a,c);
revk=revk*dot(a,c);




end