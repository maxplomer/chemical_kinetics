clc;clear;



A(1)=1;
B(1)=0;
C(1)=0;
k1=1;k2=1;

J=[-k1  0  0;
    k1  -k2  0
    0  k2  0];
[V,~] = eig(J)
D=eig(J)


dt=0.001;
time=0:dt:5;
for i=2:length(time)
    dAdt=-k1*A(i-1);
    dBdt=k1*A(i-1)-k2*B(i-1);
    dCdt=k2*B(i-1);
    A(i)=A(i-1)+dAdt*dt;
    B(i)=B(i-1)+dBdt*dt;
    C(i)=C(i-1)+dCdt*dt;
end
plot(time,A,time,B,time,C)

