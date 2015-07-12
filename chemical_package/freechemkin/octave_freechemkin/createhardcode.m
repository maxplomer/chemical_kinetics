function createhardcode
global fid elements species;
global composition thermdata commonT;
global nur nuf M A B E Meff Low Troe Rtype Duplicate;
addpath (fullfile(pwd,'modules'));
fid = fopen('chem.inp', 'r');
while 1
    tline = fgetl(fid);

    if(tline==-1)                   %end of file
        break;
    end
    tline=strread(tline,'%s','delimiter',' ');
    if (isempty(tline))             %deal with blank lines
       continue; 
    end
    if(strcmp(tline{1}(1),'!'))     %deal with commented lines
        continue;
    end
    switch tline{1}
        case 'ELEMENTS'
            disp('elements is here')%need to have blank/commented lines and use the end thing in all functions
            readelements;
            elements
        case 'SPECIES'
            disp('species is here')
            readspecies;
            species
        case 'THERMO'
            disp('thermo is here')
            readthermo;
            composition
            commonT
            thermdata
        case 'REACTIONS'
            disp('reactions is here')
            readreactions;
            nur
            nuf
            M
            A
            B
            E
            Meff
            Low
            Troe
            Rtype
            Duplicate
    end
end
fclose(fid);
%%%%%%%%%         Done reading chem.inp                    %%%%%%%%%
%%%%%%%%%         Now start creating hardcode files        %%%%%%%%%
if(~exist('hardcode','dir'))
    copyfile('allmech','hardcode');
end

printgetelements;
printgetcomposition;
printgetII;
printgetKK;
printgetspecies;
printgetwt;
printgetthermodata;
printgetabe;
printgetnu;
printgetnunet;
printkfkr;


end