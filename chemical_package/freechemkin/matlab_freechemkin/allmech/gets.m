function s=gets(T)
%entropy at standard state, s=f(P,T)
%units: ergs/(mole*K)      erg=g*cm^2/s^2=1e-7 J
global N commonT;
RU=83145100; %erg/(mol*K)

T2=T*T;T3=T2*T;T4=T3*T;logT=log(T);
s_over_R=zeros(getKK,1);
for i=1:length(commonT)
    if T>commonT(i)
        s_over_R(i)=N(1,i)*logT +N(2,i)*T +N(3,i)*(T2)/2 +N(4,i)*(T3)/3 +N(5,i)*(T4)/4 +N(7,i);
    else
        s_over_R(i)=N(8,i)*logT +N(9,i)*T +N(10,i)*(T2)/2 +N(11,i)*(T3)/3 +N(12,i)*(T4)/4 +N(14,i);
    end
end

s=s_over_R*RU;

end