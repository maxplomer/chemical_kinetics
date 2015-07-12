function wdot=getwc(T,Y)
%units: mol/(cm^3*sec)
%Y is mol concentration [mol/cm^3]



[fwdk,revk]=getkfkr(T,Y);
nu=getnunet; 

netk=fwdk-revk;

netk=netk';
wdot=nu*netk;

% for i=1:KK
% wdot(i)=sum(nu(i,:).*netk);
% end

%number of reactions
% II=length(netk);
% for i=1:II
%     netk=fwdk(i)-revk(i);
%     nu=nur(:,i)-nuf(:,i);
%     wdot=wdot+netk*nu;
% end

end