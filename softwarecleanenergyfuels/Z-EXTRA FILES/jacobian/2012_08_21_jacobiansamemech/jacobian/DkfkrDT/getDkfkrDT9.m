function [DfwdkDT,DrevkDT]=getDkfkrDT9(T,C)
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
logF=log(Fcent)*(1+(  (log(Pr)+c)/(n-d*(log(Pr)+c))  )^2)^-1;
F=exp(logF);


%C is molar concentration c=X*P/(RU*T)
[fwdk,revk]=ckkfkrc(T, C, chem);

Dk0DT=k0*(E0/T^2+B0/T);
DkinfDT=kinf*(E(I)/T^2 +   B(I)/T );
D1_kinfDT = -kinf^-1*(E(I)/T^2 +   B(I)/T );
DPrDT=Cm*(Dk0DT/kinf+k0*D1_kinfDT);
DPr_1_PrDT=(1+Pr)^-2*DPrDT;

XP=(log10(Pr)+c)/(n-d*(log10(Pr)+c));

DFcentDT=(1-aa)*(-1/T3)*exp(-T/T3)+aa*(-1/T1)*exp(-T/T1);
DlnFcentDT = DFcentDT/Fcent;

DnDT = -1.27*DlnFcentDT;
DcDT = -0.67*DlnFcentDT;

Dtemp1DT = DPrDT/(Pr*log(10)) + DcDT;
DXPDT =  (-(log(Pr)+c)*DnDT   +     n*Dtemp1DT)*(n-d*(log10(Pr)+c))^-2;
DlnFDT =   DlnFcentDT/(1+XP^2)     -(log(Fcent)*2*XP)/(1+XP^2)^2*DXPDT;

DfwdkDT = fwdk(I)*(DkinfDT/kinf + DPr_1_PrDT*(1+Pr)/Pr + DlnFDT);


%Kp(T)=exp{-[sum i=1:N of (vi''-vi')*mu0,i(T)]/(RU*T)}
g = ckgml(T, chem);
Kp=exp(   -dot(nu(:,I),g)  /  (RU*T)  );
%kf/kb=Kc=Kp/(RU*T)^[sum i=1:N of vi'' - vi']
Kc=Kp/(RU*T/PA)^sum(nu(:,I));

DgDT=getDgDT(T);
DKpDT=Kp*(dot(nu(:,I),g)/(RU*T^2)     - (RU*T)^-1*dot(nu(:,I),DgDT)  );
D1KcDT=-Kc^-2*(DKpDT/(RU*T/PA)^sum(nu(:,I))+Kp*(-sum(nu(:,I)))*T^(-sum(nu(:,I))-1)/(RU/PA)^sum(nu(:,I)));

%C is molar concentration c=X*P/(RU*T)


DrevkDT = revk(I)*( DkinfDT/kinf + DPr_1_PrDT*(1+Pr)/Pr + DlnFDT  +  Kc*D1KcDT );
end