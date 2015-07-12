function [DfwdkDC,DrevkDC]=getDkfkrDC18(T,C,j)
%reaction number
I=18;
chem=ckinit;
[RU, ~, PA] = ckrp(chem);
[A, B, E] = ckabe(chem); %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
                
nu=cknu(chem);
nuf=abs(cknuf(chem));
nur=nu+nuf;

kf=A(I)*T^B(I)*exp(-E(I)/T);

%C is molar concentration c=X*P/(RU*T)
DfwdkDC=kf;
for i=1:9
   if i==j
       DfwdkDC=DfwdkDC*nuf(i,I)*C(i)^(nuf(i,I)-1);
   else
       DfwdkDC=DfwdkDC*C(i)^nuf(i,I);
   end
end

%Kp(T)=exp{-[sum i=1:N of (vi''-vi')*mu0,i(T)]/(RU*T)}
g = ckgml(T, chem);
Kp=exp(   -dot(nu(:,I),g)  /  (RU*T)  );
%kf/kb=Kc=Kp/(RU*T)^[sum i=1:N of vi'' - vi']
Kc=Kp/(RU*T/PA)^sum(nu(:,I));
kb=kf/Kc;

%C is molar concentration c=X*P/(RU*T)
DrevkDC=kb;
for i=1:9
    if i==j
        DrevkDC=DrevkDC*nur(i,I)*C(i)^(nur(i,I)-1);
    else
        DrevkDC=DrevkDC*C(i)^nur(i,I);
    end
end

end