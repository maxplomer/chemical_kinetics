function [nuf,nur]=getnu
% SPECIES: H2 O2 H2O H O OH HO2 H2O2 N2 
nuf=[0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     1     1     1     1     0     0     0     0     0     1     1     0     0     0
     0     1     0     1     0     2     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     1
     0     0     0     0     0     0     0     0     0     1     1     1     1     2     2     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     1     1     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0];
 
nur=[0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     1     1     1     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0     0     1     1
     0     1     1     0     2     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     2     0     0     1     0     0     0     2     1     0     0     0     2     1     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     1     1     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0];

end