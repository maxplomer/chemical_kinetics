function createhardcode(inputfiles)
%valid input file lists are: {'chem.inp','term.dat'} or {'chem.inp'}
global thermdatexists fid fid2 elements species;
global composition thermdata commonT;
global nur nuf M A B E Meff Low Troe Rtype Duplicate;
addpath (fullfile(pwd,'modules'));
switch length(inputfiles)
    case 1
        thermdatexists=0;
        disp('therm.dat does not exist');
        fid = fopen(inputfiles{1}, 'r');
    case 2
        thermdatexists=1;
        disp('therm.dat does exist');
        fid = fopen(inputfiles{1}, 'r');
        fid2 = fopen(inputfiles{2}, 'r');
    otherwise
        disp('not valid input file list');
        return;
end
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

if (thermdatexists==1)
    disp('now going to call readthermo.m');%it will tline read from fid or fid2 depending on global thermdatexists value, also need its output to be global variables
    readthermo;
    composition
    commonT
    thermdata
end

fclose(fid);
if (thermdatexists==1)
    fclose(fid2);
end
%%%%%%%%%Done reading chem.inp (and therm.dat if specified)%%%%%%%%%
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
printkfkr;

end