function readthermo
global thermdatexists fid fid2 elements species;
global composition thermdata commonT;

thermdata=[];
composition=zeros(length(species),length(elements));
phase=[];
lowT=zeros(length(species),1);
highT=zeros(length(species),1);
commonT=[];
e=[]; %element in composition
n=[]; %element coefficient number

while 1
    if (thermdatexists==0)
        tline = fgetl(fid);
    else
        tline = fgetl(fid2);
    end
    
    name=strread(tline,'%s','delimiter',' ');
    if (isempty(name))             %deal with blank lines
       continue; 
    end
    if(strcmp(name{1}(1),'!'))     %deal with commented lines
        continue;
    end
    if (strcmp(name{1},'END'))     %deal with end keyword
        break;
    end
    name=name{1};%get species name, in first line (of 4) info
    
    
    I_species=find(strcmp(species, name));%will automatically identify the next 3 lines as "not in mech lines" if necessary
    if(isempty(I_species))
        continue;
    end
    
    comp=tline(25:44);%sometimes put # in error in 44th column, but from what i saw it was only a zero used
    e=[];
    n=zeros(4,1); %new chemkin can have 5 elements, fix later, they put the extra one after the "1" at the end
    while 1 %tricky algorithm to deal with 1-4 elements being present
        e{1}=strread(comp(1:2),'%s','delimiter',' ');
        n(1)=str2num(comp(4:5));
        e{2}=strread(comp(6:7),'%s','delimiter',' ');
        if (isempty(e{2})) 
            break;
        end
        n(2)=str2num(comp(9:10));
        e{3}=strread(comp(11:12),'%s','delimiter',' ');
        if (isempty(e{3})) 
            break;
        end
        n(3)=str2num(comp(14:15));
        e{4}=strread(comp(16:17),'%s','delimiter',' ');
        if (isempty(e{4})) 
            break;
        end
        n(4)=str2num(comp(19:20));
        break;
    end
    for i=1:4
        if (isempty(e{i})) 
            break;
        end
        I_element=find(strcmp(elements, e{i}));
        if(~isempty(I_element))
            composition(I_species,I_element)=n(i);
        else
            disp('error could not find element');
            return;
        end
    end
    phase{I_species}=tline(45);
    lowT(I_species)=str2num(tline(46:55));
    highT(I_species)=str2num(tline(56:65));
    commonT{I_species}=tline(66:73);
    %get next 3 line info
    tline = fgetl(fid);
    thermdata{I_species,1}=tline(1:15);
    thermdata{I_species,2}=tline(16:30);
    thermdata{I_species,3}=tline(31:45);
    thermdata{I_species,4}=tline(46:60);
    thermdata{I_species,5}=tline(61:75);
    tline = fgetl(fid);
    thermdata{I_species,6}=tline(1:15);
    thermdata{I_species,7}=tline(16:30);
    thermdata{I_species,8}=tline(31:45);
    thermdata{I_species,9}=tline(46:60);
    thermdata{I_species,10}=tline(61:75);
    tline = fgetl(fid);
    thermdata{I_species,11}=tline(1:15);
    thermdata{I_species,12}=tline(16:30);
    thermdata{I_species,13}=tline(31:45);
    thermdata{I_species,14}=tline(46:60);
end

end

