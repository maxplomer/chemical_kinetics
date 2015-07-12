function [fwdk, revk] = ckkfkr(P, T, X, chem)
%
% function [fwdk, revk] = ckkfkr(P, T, X, chem)
%  Purpose:  calculate forward/reverse rate of reactions given P, T, X 
%  Input:  
%        P:   pressure
%        T:   temperature
%        X:   mole fraction 
%        chem: chemkin workspace 
%  Output: 
%        fwdk:   forward rate 
%        revk:   reverse rate 
