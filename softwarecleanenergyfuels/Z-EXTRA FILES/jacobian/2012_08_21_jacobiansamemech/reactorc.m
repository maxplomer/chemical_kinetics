function reactorc
clc;clear;format long e;
chem=ckinit;
%Ru is erg/(mol*K)
%PA is 1 atm = 1 013 250 erg / (centimeter^3)
[RU, ~, PA] = ckrp(chem);
W=ckwt(chem);
% number of variables
KK=chem.KK;
NV=KK+1;

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 800;
Phi = 1;

%variables
nC = 1:KK;
nT = KK+1;

I_H2 = ckfindsp('H2',chem);I_H2=I_H2(1);
I_O2 = ckfindsp('O2',chem);I_O2=I_O2(1);
I_N2 = ckfindsp('N2',chem);I_N2=I_N2(1);

% setup the initial conditions
X0 = zeros(KK,1);

%initial mol fractions
X0(I_H2) = Phi*2;
X0(I_O2) = 1;
X0(I_N2) = 3.76;
X0 = X0/sum(X0);



% IC vectors  into the main IC vector
c0 = zeros(NV,1);
c0(nC) = X0*P0/(RU*T0); %mol frac ->c_i
c0(nT) = T0; %temperature

%time span seconds
tspan = [0 30]; 

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-5, 'AbsTol', 1.0E-9);
[t ,c] = ode15s(@ignfun, tspan, c0, options);



% plot the results
 figure

semilogx(t, c(:, nT))
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim(tspan);



%print 2nd point to screen
% c(2, nC)'
% c(2, nT)

%print jacobian and eigenvalues to screen for 2nd point - fortran
% J=ajactvfixed(c(2, nT),c(2, nC)',chem);
% %J
% for i=1:10
%     for j=1:10
%     	J(i,j)=spa_sf(J(i,j),5); %truncate to 5 sig figs
%     end
% end
% eig(J)


%print jacobian and eigenvalues to screen for 2nd point - matlab
% J=getjacobian(c(2, nT),c(2, nC)',chem);
% %J
% for i=1:10
%     for j=1:10
%     	J(i,j)=spa_sf(J(i,j),5); %truncate to 5 sig figs
%     end
% end
% eig(J)




%save max-eigenvalue to .mat file
for k=1:length(c(:, nT))
    J=ajactvfixed(c(k, nT),c(k, nC)',chem);
    for i=1:10
        for j=1:10
        	J(i,j)=spa_sf(J(i,j),5); %truncate to 5 sig figs
        end
    end
    maxeig(k)=max(eig(J));
end

save


    % RHS functions of the ODEs
    %   input:  t is the time 
    %           c: the vector of indepdent variables 
    %   output: dydt: the vector of LHS of the ODE
    %Constant Volume and U, always have to fix 2 parameters
    function dcdt = ignfun(t, c)%time, mol concentration
        dcdt = zeros(NV,1);

        C = c(nC);
        T = c(nT);
       
        cv=ckcvml(T,chem);%vector of cv mol based
        u=ckuml(T,chem);%energy mol specific

        wdot=ckwc(T, C,chem);%molelar production rate
        % change in mol concentration per time, dC/dt
        dCdt = wdot;
           
        % change in temperature per time, dT/dt
        dTdt  = -dot(u, wdot)/dot(C, cv);
        dcdt(nC) = dCdt;
        dcdt(nT) = dTdt;
    end




end