
function cemah2
clc;clear;

%IC
y0(1)=0.5;
y0(2)=0.5;
y0(3)=0;
y0(4)=0;
y0(5)=0;
y0(6)=0;
y0(7)=0;
y0(8)=0;
y0(9)=0;
%time span 0 to 10 seconds
tspan = [0 10]; %

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);

% plot the results
figure
plot(t,y(:,1),'.',t,y(:,2),'.',t,y(:,3),'.');
xlabel('Time (sec)');
ylabel('[y]');
xlim([0 10]);

%
% for i=1:length(t)
%     J=[-16*k1*y(i,1)^3        0      0;
%          4*k1*y(i,1)^3       -k2     0;
%           0                 k2     0];
%       eigs(i,1:3)=eig(J);
% end
% figure
% plot(t,eigs(:,1),t,eigs(:,2),t,eigs(:,3))
    function dydt = ignfun(t, y)
        dydt(1)=-4*k1*y(1)^4;
        dydt(2)=k1*y(1)^4-k2*y(2);
        dydt(3)=k2*y(2);
        dydt=dydt';
    end


end


