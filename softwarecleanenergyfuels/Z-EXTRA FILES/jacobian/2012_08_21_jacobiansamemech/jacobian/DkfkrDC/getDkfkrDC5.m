function [DfwdkDC,DrevkDC]=getDkfkrDC5(T,C,j)
%reaction number
I=5;
chem=ckinit;
[RU, ~, PA] = ckrp(chem);
[A, B, E] = ckabe(chem); %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
                
nu=cknu(chem);
nuf=abs(cknuf(chem));
nur=nu+nuf;

kf=A(I)*T^B(I)*exp(-E(I)/T);

%third body efficiences
%H2/2.5/ H2O/12/
a=ones(9,1);a(1)=2.5;a(5)=12;
Cm=dot(a,C);

%C is molar concentration c=X*P/(RU*T)
DfwdkDC=kf;
for i=1:9
   if i==j
       DfwdkDC=DfwdkDC*nuf(i,I)*C(i)^(nuf(i,I)-1);
   else
       DfwdkDC=DfwdkDC*C(i)^nuf(i,I);
   end
end
[fwdk,revk]=ckkfkrc(T, C, chem);
DfwdkDC=DfwdkDC*Cm+a(j)*fwdk(I)/Cm;

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
DrevkDC=DrevkDC*Cm+a(j)*revk(I)/Cm;
end