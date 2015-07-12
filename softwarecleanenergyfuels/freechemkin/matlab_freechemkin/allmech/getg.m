function g=getg(T)
%gibbs free energy at standard state(p=1atm)
%units: erg/mol
%g=h - Ts
global N commonT;
RU=83145100; %erg/(mol*K)
%g=geth(T) - T*gets(T);

T2=T*T;T3=T2*T;T4=T3*T;T5=T4*T;
logT=T*(1-log(T));
g_over_R=zeros(getKK,1);
for i=1:length(commonT)
    if T>commonT(i)
        g_over_R(i)=N(1,i)*logT +N(2,i)*-(T2)/2 +N(3,i)*-(T3)/6 +N(4,i)*-(T4)/12 +N(5,i)*-(T5)/20 +N(6,i)-N(7,i)*T;
    else
        g_over_R(i)=N(8,i)*logT +N(9,i)*-(T2)/2 +N(10,i)*-(T3)/6 +N(11,i)*-(T4)/12 +N(12,i)*-(T5)/20 +N(13,i)-N(14,i)*T;
    end
end

g=g_over_R*RU;



end