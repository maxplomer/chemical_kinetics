function run_reactor_perfectlystirred
if exist('reactor_perfectlystirred','file')==2
    delete('reactor_perfectlystirred')	
end
if exist('output.txt','file')==2
    delete('output.txt')	
end

system('f95 -o reactor_perfectlystirred reactor_perfectlystirred.f mechspecific.f');
system('./reactor_perfectlystirred>>output.txt');

if exist('reactor_perfectlystirred','file')==2
    delete('reactor_perfectlystirred')	
end


end
