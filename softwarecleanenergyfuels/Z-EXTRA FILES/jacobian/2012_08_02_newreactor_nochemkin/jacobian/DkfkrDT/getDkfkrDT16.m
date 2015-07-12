function [DfwdkDT,DrevkDT]=getDkfkrDT16(T,C)
%reaction number

I=16;


RU=83145100; %erg/(mol*K)
[A,B,E]=getabe; %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
[nuf,nur]=getnu; 
%third body efficiences
%H2/2.5/ H2O/12/
a=ones(9,1);a(1)=2.5;a(3)=12;
Cm=dot(a,C);
%lindemann approach
A0=1.202E+17;B0=0;E0=4.55E+04;
k0=(A0)*T^(B0)*exp(-(E0)*41840000/(RU*T));
kinf=A(I)*T^B(I)*exp(-E(I)*41840000/(RU*T));
Pr=(k0*Cm)/kinf;

%Troe form
a=0.8;%alpha
T3=1E-30;%T***
T1=1E+30;%T*
Fcent=(1-a)*exp(-T/T3)+a*exp(-T/T1);   %T** not included
c = -0.4 - 0.67*log10(Fcent);
n = 0.75 - 1.27*log10(Fcent);
d = 0.14;
logF=log(Fcent)*(1+(  (log(Pr)+c)/(n-d*(log(Pr)+c))  )^2)^-1;
F=exp(logF);


%C is molar concentration c=X*P/(RU*T)
[fwdk,revk]=getkfkr(T,C);

Dk0DT=k0*(E0*41840000/(RU*T^2)+B0/T);
DkinfDT=kinf*(E(I)*41840000/(RU*T^2) +   B(I)/T );
D1_kinfDT = -kinf^-1*(E(I)*41840000/(RU*T^2) +   B(I)/T );
DPrDT=Cm*(Dk0DT/kinf+k0*D1_kinfDT);
DPr_1_PrDT=(1+Pr)^-2*DPrDT;

XP=(log10(Pr)+c)/(n-d*(log10(Pr)+c));

DFcentDT=(1-a)*(-1/T3)*exp(-T/T3)+a*(-1/T1)*exp(-T/T1);
DlnFcentDT = DFcentDT/Fcent;

DnDT = -1.27*DlnFcentDT;
DcDT = -0.67*DlnFcentDT;

Dtemp1DT = DPrDT/(Pr*log(10)) + DcDT;
DXPDT =  (-(log(Pr)+c)*DnDT   +     n*Dtemp1DT)*(n-d*(log10(Pr)+c))^-2;
DlnFDT =   DlnFcentDT/(1+XP^2)     -(log(Fcent)*2*XP)/(1+XP^2)^2*DXPDT;

DfwdkDT = fwdk(I)*(DkinfDT/kinf + DPr_1_PrDT*(1+Pr)/Pr + DlnFDT);


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


%DrevkDT=F*(Pr/(1+Pr))*DrevkDT + F*revk(I)*DPr_1_PrDT + (Pr/(1+Pr))*revk(I)*DFDT;

DrevkDT = revk(I)*( DkinfDT/kinf + DPr_1_PrDT*(1+Pr)/Pr + DlnFDT  +  Kc*D1KcDT );
end