function abcreactor
clc;clear;
%simple model

%A-k1->B
%B-k2->C

%dA/dt=-k1*A
%dB/dt=k1*A-k2*B
%dC/dt=k2*B

%reaction rates
k1=1;
k2=1;

%intial concentrations
A=1;
B=0;
C=0;

y0= [A,B,C];

%time span seconds
tspan = [0 30]; 

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-5, 'AbsTol', 1.0E-9);
[t ,y] = ode23s(@ignfun, tspan, y0, options);

% plot the results
 figure

plot(t, y(:, 1))
hold on;
plot(t, y(:, 2),'.')
hold on;
plot(t, y(:, 3),'x')


	function dydt = ignfun(t, y)
        dydt = zeros(3,1);%do i need this line? or is it just making everything zero?
k1=1;
k2=1;
        dydt(1)=-k1*y(1);
        dydt(2)=k1*y(1)-k2*y(2); %should set up with wdot function, and nu function
        dydt(3)=k2*y(2);
    end

end
