function [DWDC]=getDWDC(T,C)
[nuf,nur]=getnu;
nu=nur-nuf;
DWDC=zeros(9,9);
for i=1:9
    for j=1:9
        [DfwdkDC,DrevkDC]=getDkfkrDC(T,C,j);
        DnetkDC=DfwdkDC-DrevkDC;
        DWDC(i,j)=dot(nu(i,:),DnetkDC);         % derivative of g_i wrt species j
    end
end

end