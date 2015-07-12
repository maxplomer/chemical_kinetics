function [wdot] = ckwc(T, C, chem)
%
% function [wdot] = ckwc(T, C, chem)
%  Purpose:  calculate mole production rate of species given T, C 
%  Input:  
%        T:   temperature
%        C:   mole concentration 
%        chem: chemkin workspace 
%  Output: 
%        wdot:   mole production rate 
%
