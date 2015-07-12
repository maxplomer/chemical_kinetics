

function csppaperstudyplot_ode
clc;clear;
y0(1)=1;
y0(2)=0;
y0(3)=0;
k1=1;k2=1;

J=[-k1  0  0;
    k1  -k2  0
    0  k2  0];
[V,~] = eig(J)
d = eig(J)


%time span 0 to 10 seconds
tspan = [0 10]; %

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);

% plot the results
% figure
% plot(t,y(:,1),'.',t,y(:,2),'.',t,y(:,3),'.');
% xlabel('Time (sec)');
% ylabel('[y]');
% xlim([0 10]);


    function dydt = ignfun(t, y)
        dydt(1)=-k1*y(1);
        dydt(2)=k1*y(1)-k2*y(2);
        dydt(3)=k2*y(2);
        dydt=dydt';
    end
end