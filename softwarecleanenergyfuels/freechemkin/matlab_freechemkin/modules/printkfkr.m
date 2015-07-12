function printkfkr
global nuf nur M Meff Low Troe Rtype species A B E;
RU=83145100; %erg/(mol*K)
nu=nur-nuf;
[~,n] = size(nuf);
fid = fopen('hardcode/getkfkr.m', 'wt');
fprintf(fid, 'function [fwdk,revk]=getkfkr(T,Y)\n');
fprintf(fid, 'RU=83145100; %%erg/(mol*K)\n');
fprintf(fid, 'PATM=1.01325D6/RU;\n');
fprintf(fid, 'g=getg(T)/RU;\n');
fprintf(fid, 'sumY=sum(Y);\n');
fprintf(fid, 'revk=zeros(1,%g);\n',n);
fprintf(fid, 'fwdk=zeros(1,%g);\n',n);
fprintf(fid, '[A, B, E]=getabe;\n');
fprintf(fid, 'kfarray=A.*T.^B.*exp(-E*41840000/(RU*T));\n\n');

for i=1:n
    %put this into its own function or series of if statement
    if (isempty(M{i}) || strcmpi(M{i},'M'))%if no pressure dependence, this is for normal arrhenius or 3rd body
        fprintf(fid, 'kf=kfarray(%g)',i);
        if (strcmpi(M{i},'M'))
            fprintf(fid, '*(sumY');
            for j=1:2:(length(Meff{i})-1)
                I_species=find(strcmp(species,Meff{i}{j}));
                fprintf(fid, '+Y(%g)',I_species);
                fprintf(fid, '*%g',str2num(Meff{i}{j+1})-1);
            end
            fprintf(fid, ')');
        end
        fprintf(fid, ';\n');
        
        fprintf(fid, 'fwdk(%g)=kf',i);
        for j=1:length(species)
            if (nuf(j,i)>0)
                for k=1:nuf(j,i)
                    fprintf(fid, '*Y(%g)',j);
                end
            end
        end
        fprintf(fid, ';\n');
        
        if (strcmp(Rtype{i},'norm'))
            fprintf(fid, 'Kc=exp(   -(');
            for j=1:length(species)
                if (nu(j,i)~=0);
                    fprintf(fid, '+g(%g)*%g',j,nu(j,i));
                end
            end
            fprintf(fid, ')  /  T  )');
            if (sum(nu(:,i))~=0)
                if (sum(nu(:,i))>0)
                    for j=1:sum(nu(:,i))
                        fprintf(fid, '/(T/PATM)');
                    end
                else
                    for j=1:abs(sum(nu(:,i)))
                        fprintf(fid, '*(T/PATM)');
                    end
                end
            end
            fprintf(fid, ';\n');

            fprintf(fid, 'revk(%g)=kf/Kc',i);
            for j=1:length(species)
                if (nur(j,i)>0)
                    for k=1:nur(j,i)
                        fprintf(fid, '*Y(%g)',j);
                    end
                end
            end
            fprintf(fid, ';\n');
        end
        
        fprintf(fid, '\n');
        continue;
    end
    %%%%%troe Pressuredependence%%%%%start

    fprintf(fid, '%%lindemann approach\n');
    abeLow=strread(Low{i},'%s','delimiter',' ');
    
    
    fprintf(fid, 'k0=%s',abeLow{1});
    if (str2num(abeLow{2})~=0)
        if (str2num(abeLow{2})==1)
            fprintf(fid, '*T');
        else
            fprintf(fid, '*T^%s',abeLow{2});
        end
    end
    if (str2num(abeLow{3})~=0)
        fprintf(fid, '*exp(%e/T)',-str2num(abeLow{3})*41840000/RU);
    end
    fprintf(fid, ';\n');
    fprintf(fid, 'kinf=kfarray(%g);\n',i);
    
    fprintf(fid, 'Pr=(k0*(sumY');
    for j=1:2:(length(Meff{i})-1)
        I_species=find(strcmp(species,Meff{i}{j}));
        fprintf(fid, '+Y(%g)',I_species);
        fprintf(fid, '*%g',str2num(Meff{i}{j+1})-1);
    end
    fprintf(fid, '))/kinf;\n');

    %%%if after delimited troeform legnth =3 then do this if 4 then do this
    if (~isempty(Troe{i}))
        fprintf(fid, '%%Troe form\n');
        aTTT=strread(Troe{i},'%s','delimiter',' '); 
        fprintf(fid, 'a=%s;%%alpha\n',aTTT{1});
        fprintf(fid, 'T3=%s;%%T***\n',aTTT{2});
        fprintf(fid, 'T1=%s;%%T*\n',aTTT{3});
        switch length(aTTT)
            case 3
            	fprintf(fid, 'Fcent=(1-a)*exp(-T/T3)+a*exp(-T/T1);   %%T** not included\n');
            case 4
            	fprintf(fid, 'T2=%s;%%T**\n',aTTT{4});
            	fprintf(fid, 'Fcent=(1-a)*exp(-T/T3)+a*exp(-T/T1)+exp(-T2/T);\n');
            otherwise
                disp('not valid troe input');
                return;
        end
        fprintf(fid, 'c = -0.4 - 0.67*log10(Fcent);\n');
        fprintf(fid, 'n = 0.75 - 1.27*log10(Fcent);\n');
        fprintf(fid, 'd = 0.14;\n');
        fprintf(fid, 'F=10^(log10(Fcent)*(1+(  (log10(Pr)+c)/(n-d*(log10(Pr)+c))  )^2)^-1);\n');
    else
        fprintf(fid, 'F=1;\n');
    end
    fprintf(fid, 'kf=kinf*(Pr/(1+Pr))*F;\n');
    fprintf(fid, 'fwdk(%g)=kf',i);
    for j=1:length(species)
        if (nuf(j,i)>0)
            for k=1:nuf(j,i)
                fprintf(fid, '*Y(%g)',j);
            end
        end
    end
    fprintf(fid, ';\n');

    if (strcmp(Rtype{i},'norm'))
        fprintf(fid, 'Kc=exp(   -(');
        for j=1:length(species)
            if (nu(j,i)~=0);
                fprintf(fid, '+g(%g)*%g',j,nu(j,i));
            end
        end
        fprintf(fid, ')  /  T  )');
        if (sum(nu(:,i))~=0)
            if (sum(nu(:,i))>0)
                for j=1:sum(nu(:,i))
                    fprintf(fid, '/(T/PATM)');
                end
            else
                for j=1:abs(sum(nu(:,i)))
                    fprintf(fid, '*(T/PATM)');
                end
            end
        end
        fprintf(fid, ';\n');
        
        fprintf(fid, 'revk(%g)=kf/Kc',i);
        for j=1:length(species)
            if (nur(j,i)>0)
                for k=1:nur(j,i)
                    fprintf(fid, '*Y(%g)',j);
                end
            end
        end
        fprintf(fid, ';\n');
    end
    
    fprintf(fid, '\n');
    
    %%%%%troe Pressuredependence%%%%%end
   

end
fprintf(fid, 'end');
fclose(fid);

end

%things to remember:
%normal or forward reacitons, this info is inputed
%rev from princetonmech, have to add that keyword in
%pressure depenedence of species w/o troe, this info is inputed

%if troe or pressure dependence kf is kinf

%now need to look at how reg/Pdep/troe differs and how can setup this
    %code