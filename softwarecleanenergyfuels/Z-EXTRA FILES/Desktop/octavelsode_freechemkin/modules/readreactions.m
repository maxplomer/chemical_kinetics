function  readreactions
global  fid species;
global nur nuf M A B E Meff Low Troe Rtype Duplicate;


nur=[];     % (reaction,species)
nuf=[];     % nuf is reactants, nur is products
M=[];       % (reaction) 3rd body, if blank none, then either M ,(+M) or (+species)
A=[];
B=[];
E=[];
Meff=[];    % 3rd body efficiencies stored array in cell style, same with low/troe but only an array of length 1
Low=[];
Troe=[];
Rtype=[];   % reactiontype: normal 'norm' or foward 'forw'
Duplicate=[];
Rnum=0;     % reactionnumber

ns={'1','2','3','4','5','6','7','8','9','0','.'};%numbers for coefficients
while 1
    tline = fgetl(fid);
    eqn1=textscan(tline,'%s');%eqn1 is tline delimiter ' ', get first line info
    eqn1=eqn1{1};


    if (isempty(eqn1))             %deal with blank lines
       continue; 
    end
    if(strcmp(eqn1{1}(1),'!'))     %deal with commented lines
        continue;
    end
    if (strcmp(eqn1{1},'END'))     %deal with end keyword
        break;
    end

    eqn=eqn1{1};%might want to move this inside first while 1

    flageqn=0;%flag if equation found
    %might need a flag, if not eqn, then delimit by  '/'
    while 1 %if there is a = <=> <= or => then reactionnum=reaction+1;, couldn't convert this to a switch statement because using findstr function
        if(findstr(eqn, '<=>'))
            Rnum=Rnum+1;flageqn=1;
            Rtype{Rnum}='norm';M{Rnum}=[];Troe{Rnum}=[];
            eqn=strread(eqn,'%s','delimiter','<=>');
            eqnf=eqn{1};
            eqnr=eqn{4};

            [eqnf,eqnr]=lookforPdependence(eqnf,eqnr,Rnum);%look for pressure dependence, remove (+___) then add to M array
            lookforcoeffients(eqnf, eqnr,Rnum,ns,eqn1);
            break;
        end

        if(findstr(eqn, '=>'))
            Rnum=Rnum+1;flageqn=1;
            Rtype{Rnum}='forw';M{Rnum}=[];Troe{Rnum}=[];

            eqn=strread(eqn,'%s','delimiter','=>');
            eqnf=eqn{1};
            eqnr=eqn{3};

            [eqnf,eqnr]=lookforPdependence(eqnf,eqnr,Rnum);%look for pressure dependence, remove then add to M array
            lookforcoeffients(eqnf, eqnr,Rnum,ns,eqn1);
            break;
        end
        if(findstr(eqn, '='))
            Rnum=Rnum+1;flageqn=1;
            Rtype{Rnum}='norm';M{Rnum}=[];Troe{Rnum}=[];

            eqn=strread(eqn,'%s','delimiter','=');
            eqnf=eqn{1};
            eqnr=eqn{2};

            [eqnf,eqnr]=lookforPdependence(eqnf,eqnr,Rnum);%look for pressure dependence, remove then add to M array
            lookforcoeffients(eqnf, eqnr,Rnum,ns,eqn1);
        end
        break;
    end
    if (flageqn==0)%delimit by /, then look at first term, while 1, if troe, ..... then break, do both troe and mefficiencies as cell in array, so be no meffs for pdenpence of species


        eqn2=strread(tline,'%s','delimiter','/');%eqn2 is list of meff, troe or low
        while 1
            if(findstr(eqn2{1}, 'TROE'))
                Troe{Rnum}=eqn2{2};
                break;
            end
            if(findstr(eqn2{1}, 'LOW'))
                Low{Rnum}=eqn2{2};
                break;
            end
            if(findstr(eqn2{1}, 'DUPLICATE'))
                Duplicate(Rnum)=1;
                break;
            end
            Meff{Rnum}=eqn2;
            break;
        end

    end
end

    function [eqnf,eqnr]=lookforPdependence(eqnf,eqnr,Rnum)
	global M;
        if(strfind(eqnf, '(+'))
            k=findstr(eqnf, '(+');
            M{Rnum}=eqnf(k:length(eqnf));
            eqnf=eqnf(1:(k-1));
        end
        if(strfind(eqnr, '(+'))
            k=findstr(eqnr, '(+');
            M{Rnum}=eqnr(k:length(eqnr));
            eqnr=eqnr(1:(k-1));
        end
    end

    function lookforcoeffients(eqnf, eqnr, Rnum,ns,eqn1)
	global nuf nur M species A B E;
        eqnf=strread(eqnf,'%s','delimiter','+');
        eqnr=strread(eqnr,'%s','delimiter','+');

        nuf(:,Rnum)=zeros(length(species),1);
        nur(:,Rnum)=zeros(length(species),1);
        
        %forward reaction
        for i=1:length(eqnf)
            flag=0;
            lc=0;%coefficient lenght
            for k=1:length(eqnf{i})
                for j=1:length(ns)
                    if(strcmp(eqnf{i}(k),ns{j}))
                        lc=lc+1;
                        break;
                    else
                        flag=1;
                    end
                end
                if(flag==1)
                    break;
                end
            end
            spec=eqnf{i};
            specL=length(eqnf{i});
            
            if (lc==0)
                coeff='1';
            else
                coeff=eqnf{i}(1:lc);
                spec=spec((lc+1):specL);
            end
            
            I_species=find(strcmp(spec,species));
            if(~isempty(I_species))
                nuf(I_species,Rnum)=str2num(coeff)+nuf(I_species,Rnum);
            else
                M{Rnum}=spec;
            end
        end
            %reverse reaction
        for i=1:length(eqnr)
            flag=0;
            lc=0;%coefficient lenght
            for k=1:length(eqnr{i})
                for j=1:length(ns)
                    if(strcmp(eqnr{i}(k),ns{j}))
                        lc=lc+1;
                        break;
                    else
                        flag=1;
                    end
                end
                if(flag==1)
                    break;
                end
            end
            spec=eqnr{i};
            specL=length(eqnr{i});
            
            if (lc==0)
                coeff='1';
            else
                coeff=eqnr{i}(1:lc);
                spec=spec((lc+1):specL);
            end
            
            I_species=find(strcmp(spec,species));
            if(~isempty(I_species))
                nur(I_species,Rnum)=str2num(coeff)+nur(I_species,Rnum);
            else
                M{Rnum}=spec;
            end
        end
        A{Rnum}=eqn1{2};B{Rnum}=eqn1{3};E{Rnum}=eqn1{4};
    end
end


%Random programming notes, should be deleted:
%if find troe...then break!!!!we can just store the low,troe,meff
        %strings into their respective spots, deal with later
        
        %if find low...then break
        
        %no if statment just search for species and take position after it
        %into 3rd body





%lot of empty space might want a thing where just store coordinates and
%value

%instead of while with breaks, could do a switch/case statement

%later worry about upper/lower case in M array


%need to split away coefficient
        %need to find column that letter starts and not "123etc or ."
        %or try to search for species and delimit with species name to
        %separate
        %eqnf{1}(1)
        %length(eqnf{1})
        %str2num(eqnf(1)) %note will need for loop for the length of eqnf
        %probably want to make this little section that reads the ceofs and
        %nuf nur into a function with inputs Rnum and eqn, put after
        %Rnum=Rnum+1;   Also do tline as input and assign ABE in there
        %might need to leave lines above this out of the function because
        %delimter needs {1} and {4} spot
        %need to put into function also have input being 'f' or 'r' to know
        %whether to add to foward or reverse reaction, 
        
        
        
        %will later need to go through list of species, and then list of if () or M or
        %(M)
        
        
        %!!!!we can an issue with Pdepedent H+O2(+M)=HO2(+M) , is not
        %delimited by +, will need an extra search before hand that looks
        %for '(+' in a string (and ')'), preferably that outputs the position
        
