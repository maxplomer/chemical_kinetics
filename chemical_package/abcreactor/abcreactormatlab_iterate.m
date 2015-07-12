function abcreactormatlab_iterate
%simple model

%A-k1->B
%B-k2->C

%dA/dt=-k1*A
%dB/dt=k1*A-k2*B
%dC/dt=k2*B

%reaction rates
k1=1;      k2=1;
%intial concentrations
A=1;  B=0;   C=0;


dt=0.01;
t=0;
t_plot=t;

for i=1:800
    t=t+dt;
    t_plot=[t_plot t];
    dAdt=-k1*A(i);
    dBdt=k1*A(i)-k2*B(i);
    dCdt=k2*B(i);
    
    A=[A A(i)+dAdt*dt];
    B=[B B(i)+dBdt*dt];
    C=[C C(i)+dCdt*dt];
end


% plot the results
 figure

plot(t_plot, A)
hold on;
plot(t_plot, B)
hold on;
plot(t_plot, C)

end

