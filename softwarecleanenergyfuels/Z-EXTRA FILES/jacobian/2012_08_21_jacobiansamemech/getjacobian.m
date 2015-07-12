function [JC] = getjacobian(T, C, chem)
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
addpath (fullfile(pwd,'jacobian'));
DWDT=getDWDT(T,C);
DWDC=getDWDC(T,C);  


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