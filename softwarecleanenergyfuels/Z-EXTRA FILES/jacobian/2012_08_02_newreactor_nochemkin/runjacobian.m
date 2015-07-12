%pt160
% clc;clear;
chem=ckinit;
C=[0.000000045068517
   0.000000022534258
                   1e-22
                   1e-22
                   1e-22
                   1e-22
                   1e-22
                   1e-22
   0.000000084728811];
C=C*1.0e+002;

T= 800;

format long e


% J=ajactvmax(T, C, chem);
% J(10,:)

J=getjacobian(T,C);
J





% J=ajactv(T, C, chem);
% 
% specme={'H2','O2','H2O','H','O','OH','HO2','H2O2','N2'};
% specrq={'H2','O2','O','OH','H2O','H','HO2','H2O2','N2'};
% 
% for i=1:9
%     for j=1:9
%         Jnew(i,j)=J( find(strcmp(specrq,specme{i}))  ,  find(strcmp(specrq,specme{j})) );
%     end
% end
% 
% for i=1:9
%    Jnew(i,10)= J( find(strcmp(specrq,specme{i}))  ,  10 );
% end
% 
% Jnew
% %eig(J)

% SPECIES
%  
% !jacobian: 