function [DfwdkDC,DrevkDC]=getDkfkrDC9(T,C,j)
%reaction number
I=9;
chem=ckinit;
[RU, ~, PA] = ckrp(chem);
[A, B, E] = ckabe(chem); %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
                
nu=cknu(chem);
nuf=abs(cknuf(chem));
nur=nu+nuf;  
%third body efficiences
%H2/2.5/ H2O/12/
a=ones(9,1);a(1)=2;a(5)=11;a(2)=0.78;
Cm=dot(a,C);
%lindemann approach
A0=6.366E+20;B0=-1.72;E0=5.248E+02*41840000/RU;
k0=(A0)*T^(B0)*exp(-E0/T);
kinf=A(I)*T^B(I)*exp(-E(I)/T);
Pr=(k0*Cm)/kinf;
%Troe form
aa=0.8;%alpha
T3=1E-30;%T***
T1=1E+30;%T*
Fcent=(1-aa)*exp(-T/T3)+aa*exp(-T/T1);   %T** not included
c = -0.4 - 0.67*log10(Fcent);
n = 0.75 - 1.27*log10(Fcent);
d = 0.14;
logF=log(Fcent)*(1+(  (log10(Pr)+c)/(n-d*(log10(Pr)+c))  )^2)^-1;
F=exp(logF);

%C is molar concentration c=X*P/(RU*T)
DfwdkDC=kinf;
for i=1:9
   if i==j
       DfwdkDC=DfwdkDC*nuf(i,I)*C(i)^(nuf(i,I)-1);
   else
       DfwdkDC=DfwdkDC*C(i)^nuf(i,I);
   end
end
[fwdk,revk]=ckkfkrc(T, C, chem);

XP=(log10(Pr)+c)/(n-d*(log10(Pr)+c));
Dlog10PrDc=1/(Pr*log(10))*k0/kinf*a(j);
DXPDC=(n/(n-d*(log10(Pr)+c))^2)*Dlog10PrDc;
DlnFDC = -(log(Fcent)*2*XP)/(1+XP^2)^2*DXPDC;

DfwdkDC=DfwdkDC*F*(Pr/(1+Pr)) + fwdk(I)*(a(j)/(Cm*(1+Pr)) +    DlnFDC) ;

%Kp(T)=exp{-[sum i=1:N of (vi''-vi')*mu0,i(T)]/(RU*T)}
g = ckgml(T, chem);
Kp=exp(   -dot(nu(:,I),g)  /  (RU*T)  );
%kf/kb=Kc=Kp/(RU*T)^[sum i=1:N of vi'' - vi']
Kc=Kp/(RU*T/PA)^sum(nu(:,I));
kb=kinf/Kc;

%C is molar concentration c=X*P/(RU*T)
DrevkDC=kb;
for i=1:9
    if i==j
        DrevkDC=DrevkDC*nur(i,I)*C(i)^(nur(i,I)-1);
    else
        DrevkDC=DrevkDC*C(i)^nur(i,I);
    end
end
DrevkDC=DrevkDC*F*(Pr/(1+Pr)) + revk(I)*(a(j)/(Cm*(1+Pr)) +    DlnFDC) ;


end