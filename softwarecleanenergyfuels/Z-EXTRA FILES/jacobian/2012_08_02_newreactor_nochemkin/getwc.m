function wdot=getwc(T,Y)
%units: mol/(cm^3*sec)
%Y is mol concentration [mol/cm^3]



[fwdk,revk]=getkfkr(T,Y);
[nuf,nur]=getnu; 

netk=fwdk-revk;
nu=nur-nuf;

wdot=zeros(9,1);
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