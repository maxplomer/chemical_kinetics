function [JC]=getjacobian(T,C)

addpath (fullfile(pwd,'jacobian'));

DWDT=getDWDT(T,C);
DWDC=getDWDC(T,C);  

JCTC=zeros(1,9);
cv=getcv(T);
u=getu(T);
wdot=getwc(T, C); %first 8 sigfigs match ckwc 


for j=1:9
    JCTC(j)=-dot(u, DWDC(:,j))/dot(C, cv)+dot(u, wdot)*cv(j)/dot(C, cv)^2;
end

DcvDT=getDcpDT(T);%note: DcvDT=DcpDT

D1_cvcDT=-dot(C, DcvDT)/[dot(C, cv)]^2;
for i=1:9
    DuWDT(i)=u(i)*DWDT(i)+wdot(i)*cv(i);
end
JCTT=-sum(DuWDT)/dot(C, cv)-dot(u, wdot)*D1_cvcDT;


% assemble JC
JC = [DWDC DWDT; ...
      JCTC JCTT];  
end