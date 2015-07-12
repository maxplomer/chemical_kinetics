function [DfwdkDT,DrevkDT]=getDkfkrDT2(T,C)
%reaction number
I=2;
chem=ckinit;
[RU, ~, PA] = ckrp(chem);
[A, B, E] = ckabe(chem); %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
                
nu=cknu(chem);
nuf=abs(cknuf(chem));
nur=nu+nuf;                          


%C is molar concentration c=X*P/(RU*T)
DfwdkDT=A(I)*exp(-E(I)/T)*( T^B(I)*E(I)/T^2  +  B(I)*T^(B(I)-1)  );
for i=1:9
   DfwdkDT=DfwdkDT*C(i)^nuf(i,I);
end


%Kp(T)=exp{-[sum i=1:N of (vi''-vi')*mu0,i(T)]/(RU*T)}
g = ckgml(T, chem);
Kp=exp(   -dot(nu(:,I),g)  /  (RU*T)  );
%kf/kb=Kc=Kp/(RU*T)^[sum i=1:N of vi'' - vi']
Kc=Kp/(RU*T/PA)^sum(nu(:,I));

DgDT=getDgDT(T);
DKpDT=Kp*(dot(nu(:,I),g)/(RU*T^2)     - (RU*T)^-1*dot(nu(:,I),DgDT)  );
D1KcDT=-Kc^-2*(DKpDT/(RU*T/PA)^sum(nu(:,I))+Kp*(-sum(nu(:,I)))*T^(-sum(nu(:,I))-1)/(RU/PA)^sum(nu(:,I)));


%C is molar concentration c=X*P/(RU*T)
DrevkDT=A(I)*exp(-E(I)/T)*( T^B(I)/Kc*E(I)/(T^2)  +  B(I)/Kc*T^(B(I)-1) +T^B(I)*D1KcDT );
for i=1:9
   DrevkDT=DrevkDT*C(i)^nur(i,I);
end

end