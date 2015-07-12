function [DfwdkDT,DrevkDT]=getDkfkrDT7(T,C)
%reaction number
I=7;


RU=83145100; %erg/(mol*K)
[A,B,E]=getabe; %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
[nuf,nur]=getnu;                            
%third body efficiences
%H2/2.5/ H2O/12/
a=ones(9,1);a(1)=2.5;a(3)=12;
Cm=dot(a,C);

%C is molar concentration c=X*P/(RU*T)
DfwdkDT=A(I)*exp(-E(I)*41840000/(RU*T))*( T^B(I)*(E(I)*41840000/(RU*T^2))  +  B(I)*T^(B(I)-1)  );
for i=1:9
   DfwdkDT=DfwdkDT*C(i)^nuf(i,I);
end
DfwdkDT=DfwdkDT*Cm;

%Kp(T)=exp{-[sum i=1:N of (vi''-vi')*mu0,i(T)]/(RU*T)}

g=getg(T);
Kp=exp(   -dot(nur(:,I)-nuf(:,I),g)  /  (RU*T)  );

%kf/kb=Kc=Kp/(RU*T)^[sum i=1:N of vi'' - vi']
PATM=1.01325D6;
Kc=Kp/(RU*T/PATM)^sum(nur(:,I)-nuf(:,I));

DgDT=getDgDT(T);
DKpDT=Kp*(dot(nur(:,I)-nuf(:,I),g)/(RU*T^2)     - (RU*T)^-1*dot(nur(:,I)-nuf(:,I),DgDT)  );
D1KcDT=-Kc^-2*(DKpDT/(RU*T/PATM)^sum(nur(:,I)-nuf(:,I))+Kp*(-sum(nur(:,I)-nuf(:,I)))*T^(-sum(nur(:,I)-nuf(:,I))-1)/(RU/PATM)^sum(nur(:,I)-nuf(:,I)));


%C is molar concentration c=X*P/(RU*T)
DrevkDT=A(I)*exp(-E(I)*41840000/(RU*T))*( T^B(I)/Kc*(E(I)*41840000/(RU*T^2))  +  B(I)/Kc*T^(B(I)-1) +T^B(I)*D1KcDT );
for i=1:9
   DrevkDT=DrevkDT*C(i)^nur(i,I);
end
DrevkDT=DrevkDT*Cm;
end