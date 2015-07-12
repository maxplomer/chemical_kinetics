function DgDT=getDgDT(T)
%gibbs free energy at standard state(p=1atm)
%units: erg/mol
%g=h - Ts
chem=ckinit;
[RU, ~, ~] = ckrp(chem);
DgDT=ckcpml(T, chem) - T*getDsDT(T)' - cksml(T, chem);








end