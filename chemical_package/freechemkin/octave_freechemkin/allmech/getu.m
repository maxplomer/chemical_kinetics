function u=getu(T)
%units: erg/mol
%ideal gas: p*v*M = RU*T     where v is volume/mass  and M is mass/mol
%h = u + pv    where v is volume/mol

RU=83145100; %erg/(mol*K)
u=geth(T) - RU*T;








end