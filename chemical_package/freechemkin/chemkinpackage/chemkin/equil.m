function [sol] = equil(NOP, P0, T0, TEST, X0, wrk)
%
% function [sol] = equil(NOP, P0, T0, TEST, X0, wrk)
%  Purpose:  compute chemical equilibrium 
%  Input:  
%        NOP:  problem type
%                 1: given T, P
%                 2: given T, V
%                 3: given T, S
%                 4: given P, V
%                 5: given P, H
%                 6: given P, S
%                 7: given V, U
%                 8: given V, H
%                 9: given V, S
%                10: Chapman-Jouget detonation
%        P0:   pressure (atm)
%        T0:   temperature (K)
%        TEST: guess of temperature (K)
%        X0:   inlet species mole fraction, double array of size KK
%        wrk: chemkin workspace 
%  Output: 
%        sol: solution, structure
%
