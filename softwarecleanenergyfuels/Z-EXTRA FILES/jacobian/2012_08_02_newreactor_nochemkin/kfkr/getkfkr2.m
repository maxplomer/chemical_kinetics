function [fwdk,revk]=getkfkr2(T,Y)
%reaction number
I=2;


RU=83145100; %erg/(mol*K)
[A,B,E]=getabe; %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
[nuf,nur]=getnu;                            
kf=A(I)*T^B(I)*exp(-E(I)*41840000/(RU*T));

%Y is molar concentration c=X*P/(RU*T)
fwdk=kf;
for i=1:9
   fwdk=fwdk*Y(i)^nuf(i,I);
end


%Kp(T)=exp{-[sum i=1:N of (vi''-vi')*mu0,i(T)]/(RU*T)}

g=getg(T);
Kp=exp(   -dot(nur(:,I)-nuf(:,I),g)  /  (RU*T)  );

%kf/kb=Kc=Kp/(RU*T)^[sum i=1:N of vi'' - vi']
PATM=1.01325D6;
Kc=Kp/(RU*T/PATM)^sum(nur(:,I)-nuf(:,I));
kb=kf/Kc;

%Y is molar concentration c=X*P/(RU*T)
revk=kb;
for i=1:9
   revk=revk*Y(i)^nur(i,I);
end

end