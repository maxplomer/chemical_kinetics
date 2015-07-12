%Reaction rates
k1 = 1;      k2 = 1;

%Intial concentrations
A = 1;  B=0;   C=0;

%dt is our time-step
dt = 0.01;
time = 0:dt:8;

%Iterative loop
for i=1:length(time)-1
    dAdt = -k1*A(i);
    dBdt = k1*A(i)-k2*B(i);
    dCdt = k2*B(i);
    A = [A   A(i)+dAdt*dt];
    B = [B   B(i)+dBdt*dt];
    C = [C   C(i)+dCdt*dt];
end

% plot the results
plot(time, A)
hold on;
plot(time, B ,'.')
hold on;
plot(time, C ,'x')