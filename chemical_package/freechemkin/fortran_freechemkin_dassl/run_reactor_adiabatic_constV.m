function run_reactor_adiabatic_constV
if exist('reactor_adiabatic_constV','file')==2
    delete('reactor_adiabatic_constV')	
end
if exist('output.txt','file')==2
    delete('output.txt')	
end

system('f95 -o reactor_adiabatic_constV reactor_adiabatic_constV.f mechspecific.f');
system('./reactor_adiabatic_constV>>output.txt');

if exist('reactor_adiabatic_constV','file')==2
    delete('reactor_adiabatic_constV')	
end


end
