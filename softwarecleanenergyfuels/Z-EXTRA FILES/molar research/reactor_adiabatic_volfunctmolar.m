function reactor_adiabatic_volfunctmolar

global N commonT;
[N,commonT]=getthermodata;

RU=83145100;%erg/(mol*K)
PA=1013250; %1 atm = 1 013 250 erg / (centimeter^3)
W=getwt;

% number of variables
KK=getKK;
NV=KK+2;

%system conditions
p0 = 1; % atm
P0 = p0*PA;
T0 = 500;
Phi = 0.3;
f=194.3;%rotation frequency,this is a value I will vary, original=100
R=16;%compression ratio, this is a value I will vary, original=50
omega=2*pi*f;%angular velocity
alpha=(R-1)/(R+1);

%variables
nY = 1:KK;
nT = KK+1;
nP = KK+2;

I_CH4 = findsp('CH4');
I_O2 = findsp('O2');
I_N2 = findsp('N2');

% setup the initial conditions
X0 = zeros(KK,1);

%initial mole fractions
X0(I_CH4) = Phi*0.5;
X0(I_O2) = 1;
X0(I_N2) = 3.76;
X0 = X0/sum(X0);

%initial average molecular weight
WM0=sum(X0.*W);

%initial specific volume
v0 = RU*T0 / (P0*WM0);

% here is where i put the different IC vectors 
% into the main IC vector
y0 = zeros(NV,1);
y0(nY) = X0*P0/(RU*T0);
y0(nT) = T0; %temperature
y0(nP) = P0; %pressure

%time span 0 to 10 cycle
tspan = [0 1/f]; % 1 cycles

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);

% plot the results
figure
plot(t, y(:, nT), t, y(:, 14)*20000, '.');
xlabel('Time (sec)');
ylabel('Temperature (K)');
xlim([0 1/f]);

    % RHS functions of the ODEs
    %   input:  t is the time 
    %           y: the vector of indepdent variables 
    %   output: dydt: the vector of LHS of the ODE
    %
    function dydt = ignfun(t, y)%time, mole concentration
        dydt = zeros(NV,1);

        Y = y(nY);
        T = y(nT);
        P = y(nP);

        ctot = sum(Y);
        
        cv = getcv(T);
        
        cvm = sum(cv.*Y)/ctot;
        
        u = getu(T);
        
        wdot = getwc(T, Y);
        
        %calculate change in specific volume
        X=Y/ctot;
        
        Wave = sum(W.*X);
        
        v = v0/(1+alpha)*(1+alpha*cos(omega*t));
       
        dvmassdt = v0/(1+alpha)*(-alpha*omega*sin(omega*t));
        
        dWavedt = sum(W.*wdot)/ctot;
        
        dvdt = dvmassdt*Wave+v*dWavedt;
        
        % change in mole concentration per time, dY/dt
        dYdt = wdot;
           
        % change in temperature per time, dT/dt
        dTdt = -(P*dvdt + dot(u,wdot)/ctot)/cvm;
       
        %Change in Pressure per time, dP/dt 
        dPdt = (RU*dTdt-P*dvdt)*ctot;
        
        dydt(nY) = dYdt;
        dydt(nT) = dTdt;
        dydt(nP) = dPdt;
        
    end
end
