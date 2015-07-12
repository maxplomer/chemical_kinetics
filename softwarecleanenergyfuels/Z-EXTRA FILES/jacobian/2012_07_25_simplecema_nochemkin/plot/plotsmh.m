clc;clear;
%units: erg/(mol)     erg=g*cm^2/s^2=1e-7 J
T=300:50:3000;
chem=ckinit;
[~,k_names]=ckname(chem);
for i=1:length(T)
    smh(:,i)=getsmh(T(i));
end
% for i=1:length(T)
%     cksmh(:,i)=cksmh(T(i),chem);
% end
for i=1:9
    figure
    plot(T,smh(i,:))
    hold on
%     plot(T,cksmh(i,:),'.')
    xlabel('T(K)')
    ylabel('smh')
    title(sprintf('%s',k_names{i}))
end
T=1000:50:3000;
for i=1:length(T)
    

TLOG = log(T(i));
TI = 1/T(i);

      TN(1) = TLOG - 1;
      TN(2) = T(i);
      TN(3) = TN(2)*T(i);
      TN(4) = TN(3)*T(i);
      TN(5) = TN(4)*T(i);


      SMH(i) = -1.35511017D0 +8.35033997D2*TI    ...
              +2.99142337D0*TN(1) +3.50032205D-4*TN(2)    ...
              -9.38971448D-9*TN(3) -7.69298182D-13*TN(4)   ...
              +7.91375895D-17*TN(5) 
          
end
figure
plot(T,SMH)