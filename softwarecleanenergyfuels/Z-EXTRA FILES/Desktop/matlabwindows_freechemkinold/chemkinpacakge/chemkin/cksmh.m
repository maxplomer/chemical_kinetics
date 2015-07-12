function [smh] = cksmh(T, chem)
%
% function [smh] = cksmh(T, chem)
%  Purpose:  calculate entropy of species in molar unit for give T 
%  Input:  
%        T:   temperature, K
%        chem: chemkin workspace 
%  Output: 
%        smh:   entropies minus enthalpies for the species, S(K)/R - H(K)/RT 
%
