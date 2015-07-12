function csppaperstudyplot2
clc;clear;
R=2;
k1=1;k2=1;

y(1,1)=1;
y(2,1)=0;
y(3,1)=0;



S1=[-1;1;0];
S2=[0;-1;1];



J=[-k1  0  0;
    k1  -k2  0
    0  k2  0];
[V,~] = eig(J)
D=eig(J)


dt=0.001;
time=0:dt:5;
for i=2:length(time)
    dydt=S1*F1(y(:,i-1))+S2*F2(y(:,i-1));
    y(:,i)=y(:,i-1)+ dydt*dt;
end
plot(time,y(1,:),time,y(2,:),time,y(3,:))


    function [F1]=F1(y)
        F1=k1*y(1);
    end
    function [F2]=F2(y)
        F2=k2*y(2);
    end

end
