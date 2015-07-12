function callgetkfkr
T=1000;

C=ones(1,53);
C=C';

[fwdk, revk] = getkfkr(T, C); %molelar production rate

end
