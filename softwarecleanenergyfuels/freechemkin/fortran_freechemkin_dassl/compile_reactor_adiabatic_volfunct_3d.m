function compile_reactor_adiabatic_volfunct_3d()

if exist('reactor_adiabatic_volfunct_3d','file')==2
    delete('reactor_adiabatic_volfunct_3d')	
end

system('f95 -o reactor_adiabatic_volfunct_3d reactor_adiabatic_volfunct_3d.f mechspecific.f');

end
