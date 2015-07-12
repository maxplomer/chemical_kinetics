%
% function [I] = findsp(name, splist)
%  Purpose:  find the index of a species 
%  Input:  
%        name:  name of a species, string (case sensitive)
%        splist: the list of species, StringCellArray 
%  Output: 
%        I: index of the species, empty if not found
%

function [I] = findsp(name, splist)

[trash, IA, I] = intersect({name}, splist);

end