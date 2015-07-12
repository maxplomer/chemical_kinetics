function [wdot] = ckwyp(P, T, Y, chem)
%
% function [wdot] = ckwyp(P, T, Y, chem)
%  Purpose:  calculate mole production rate of species given P, T, Y 
%  Input:  
%        P:   pressure
%        T:   temperature
%        Y:   mass fraction 
%        chem: chemkin workspace 
%  Output: 
%        wdot:   mole production rate 
%
