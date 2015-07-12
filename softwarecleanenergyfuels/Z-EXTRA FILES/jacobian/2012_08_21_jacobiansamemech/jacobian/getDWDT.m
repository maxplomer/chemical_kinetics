function [DWDT]=getDWDT(T,C)
chem=ckinit;
nu=cknu(chem);
DWDT=zeros(9,1);
for i=1:9
    [DfwdkDT,DrevkDT]=getDkfkrDT(T,C);
    DnetkDT=DfwdkDT-DrevkDT;
    DWDT(i)=dot(nu(i,:),DnetkDT);         % derivative of g_i wrt T
end

end