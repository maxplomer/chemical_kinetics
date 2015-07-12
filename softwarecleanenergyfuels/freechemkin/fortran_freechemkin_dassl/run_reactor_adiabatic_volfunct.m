function run_reactor_adiabatic_volfunct
if exist('reactor_adiabatic_volfunct','file')==2
    delete('reactor_adiabatic_volfunct')	
end
if exist('output.txt','file')==2
    delete('output.txt')	
end

system('f95 -o reactor_adiabatic_volfunct reactor_adiabatic_volfunct.f mechspecific.f');
system('./reactor_adiabatic_volfunct>>output.txt');

if exist('reactor_adiabatic_volfunct','file')==2
    delete('reactor_adiabatic_volfunct')	
end


end
