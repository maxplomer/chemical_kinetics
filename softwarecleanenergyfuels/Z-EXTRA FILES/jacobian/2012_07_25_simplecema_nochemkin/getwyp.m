function wdot=getwyp(P,T,Y)
%units: mol/(cm^3*sec)


wdot=zeros(9,1);
%initial average molecular weight
W=getwt;
WM0=1/sum(Y./W);
X=Y*WM0./W;
[fwdk,revk]=getkfkr(P,T,X);
[nuf,nur]=getnu; 


netk=fwdk-revk;
nu=nur-nuf;
for i=1:9
    wdot(i)=dot(nu(i,:),netk);
end

%number of reactions
% N=21;
% for i=1:N
%     netk=fwdk(i)-revk(i);
%     nu=nur(:,i)-nuf(:,i);
%     wdot=wdot+netk*nu;
% end

end