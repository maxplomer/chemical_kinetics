function runodesolver
if exist('odesolver','file')==2
    delete('odesolver')	
end
if exist('output.txt','file')==2
    delete('output.txt')	
end

%system('f95 -w -o odesolver odesolver.f');
system('f95 -o odesolver odesolver.f');
system('./odesolver>>output.txt');

if exist('odesolver','file')==2
    delete('odesolver')	
end


end
