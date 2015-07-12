function DgDT=getDgDT(T)
%gibbs free energy at standard state(p=1atm)
%units: erg/mol
%g=h - Ts

RU=83145100; %erg/(mol*K)
DgDT=getcp(T) - T*getDsDT(T) - gets(T);








end