function g=getg(T)
%gibbs free energy at standard state(p=1atm)
%units: erg/mol
%g=h - Ts

RU=83145100; %erg/(mol*K)
g=geth(T) - T*gets(T);








end