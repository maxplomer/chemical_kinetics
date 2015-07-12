function csppaper_Dryer_2stageign
clc;clear;
format longe;
R=3;
KK=3;
k1=1;k2=2;k3=3;

y0(1)=1;
y0(2)=0;
y0(3)=0;

S1=[-1; 1; 0];
S2=[0 ;-1; 1];
S3=[0 ;1; -1];
S=[S1  S2  S3];



J=[ -k1    0     0 ;
     k1   -k2    k3 ;
     0     k2   -k3 ];


H=[;
    ;
     ];
 
 
lamb = eig(J)





%time span 0 to 10 seconds
tspan = [0 10]; %

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);






figure
plot(t,y(:,1),t,y(:,2),'x',t,y(:,3))


    function dydt = ignfun(t, y)
        
        %dydt = S1*F1(y)+S2*F2(y)+S3*F3(y);
        
        dydt = 0;
        F=getF(y);
        for r=1:R
            dydt = dydt + S(:,r)*F(r);
        end

        



%         F=getF(y);
%         for k=1:KK
%             dydt(k) = sum(S(k,:)'.*F');
%         end
%         dydt=dydt';
    end




    function [F]=getF(y)
        F=[F1(y) F2(y) F3(y)];
    end
    function [F1]=F1(y)
        F1=k1*y(1);
    end
    function [F2]=F2(y)
        F2=k2*y(2);
    end
    function [F3]=F3(y)
        F3=k3*y(3);
    end








end
