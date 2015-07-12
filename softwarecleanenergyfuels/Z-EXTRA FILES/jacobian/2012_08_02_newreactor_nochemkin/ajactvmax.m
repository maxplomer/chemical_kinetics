function [JC] = ajactvmax(T, C, chem)
%
% function [JC] = ajactv(T, C, chem)
%  Purpose:  compute the analytic Jacobian matrix given T and C 
%
%  Input:  
%        T: temperature, K 
%        C: species mole concentration, mole/cm^3 
%        chem: chemkin workspace returned by ckinit
%  Output: 
%        jac, with temperature as the last variable
%
%  The species equations:
%       dCi_dt = FC_i = WDOT(i), i=1,KK
%  Energy equations
%       dT/dt = FT = -dot(E, WDOT)/dot(C, CV)
%

KK = chem.KK;

% compute dWDOT/d[C, T], using the analytic Jacobian
% size of J: KK by (KK+1)

specme={'H2','O2','H2O','H','O','OH','HO2','H2O2','N2'};
specrq={'H2','O2','O','OH','H2O','H','HO2','H2O2','N2'};

for i=1:9
    Cnew(i)=C(find(strcmp(specrq{i},specme)));
end

J = jact(T, Cnew);
DWDC = J(:, 1:KK);
DWDT = J(:, KK+1);
clear J

for i=1:9
    for j=1:9
        DWDCnew(i,j)=DWDC( find(strcmp(specrq,specme{i}))  ,  find(strcmp(specrq,specme{j})) );
    end
end
DWDC=DWDCnew;
for i=1:9
   DWDTnew(i)= DWDT( find(strcmp(specrq,specme{i})) );
end
DWDT=DWDTnew';
% species mole production rate, mole/cm^3/s
WDOT = ckwc(T, C, chem);

E = ckuml(T, chem); % species mole specific energy
CV = ckcvml(T, chem); % species mole specific cv
CVV = dot(C, CV); % RHO*cvbs = dot(C, CV)

% compute dFT/dC
JCTC = -E'*DWDC/CVV +dot(E, WDOT)/CVV^2*CV';

% compute dFT/dT
CVDT = ckcpdt(T, chem); % dCv/dT = dCp/dT
JCTT = -E'*DWDT/CVV -CV'*WDOT/CVV ...
       +E'*WDOT/CVV^2 * (C'*CVDT);

% assemble JC
JC = [DWDC DWDT; ...
      JCTC JCTT];  

end