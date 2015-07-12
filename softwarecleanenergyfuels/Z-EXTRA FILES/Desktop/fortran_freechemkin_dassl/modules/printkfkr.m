function printkfkr
global nuf nur M Meff Low Troe Rtype species A B E;
RU=83145100; %erg/(mol*K)
nu=nur-nuf;
[~,n] = size(nuf);
fid = fopen('odesolver.f', 'a');
fprintf(fid, '\n      SUBROUTINE GETKFKR (KK,II,T,Y,FWDK,REVK)\n');
fprintf(fid, '      INTEGER KK II\n');
fprintf(fid, '      DOUBLE PRECISION RU, PATM, T, d, sumY, kf, Kc, Z, X, k0, kinf\n');
fprintf(fid, '      DOUBLE PRECISION Pr, F, Fcent, a, T1, T2, T3, c, n\n');
fprintf(fid, '      DOUBLE PRECISION, dimension(KK) :: G, Y\n');
fprintf(fid, '      DOUBLE PRECISION, dimension(II) ::  KFARRAY, FWDK, REVK\n');
fprintf(fid, '      CALL GETRU (RU)\n');
fprintf(fid, '      CALL GETPATM (PATM)\n');
fprintf(fid, '      CALL GETG(KK,T,G)\n');
fprintf(fid, '      CALL GETKFARRAY (II, T, KFARRAY)\n');
fprintf(fid, '      PATM=PATM/RU\n');
fprintf(fid, '      G=G/RU\n');
fprintf(fid, '      sumY=SUM(Y)\n');
fprintf(fid, '      d = 0.14\n\n');


for i=1:n
    %put this into its own function or series of if statement
    if (isempty(M{i}) || strcmpi(M{i},'M'))%if no pressure dependence, this is for normal arrhenius or 3rd body
        if (strcmpi(M{i},'M'))
            fprintf(fid, '      Z=sumY\n');
            for j=1:2:(length(Meff{i})-1)
                I_species=find(strcmp(species,Meff{i}{j}));
                fprintf(fid, '      Z=Z+Y(%g)*(%g)\n',I_species,str2num(Meff{i}{j+1})-1);
            end
        end
        fprintf(fid, '      kf=KFARRAY(%g)',i);
        if (strcmpi(M{i},'M'))
           fprintf(fid, '*Z');
        end
        fprintf(fid, '\n');
        
        fprintf(fid, '      FWDK(%g)=kf',i);
        for j=1:length(species)
            if (nuf(j,i)>0)
                for k=1:nuf(j,i)
                    fprintf(fid, '*Y(%g)',j);
                end
            end
        end
        fprintf(fid, '\n');
        
        if (strcmp(Rtype{i},'norm'))
            fprintf(fid, '      X=0\n');
            for j=1:length(species)
                if (nu(j,i)~=0);
                    fprintf(fid, '      X=X+g(%g)*(%g)\n',j,nu(j,i));
                end
            end
            fprintf(fid, '      Kc=EXP(-X/T)');
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
            fprintf(fid, '\n');

            fprintf(fid, '      REVK(%g)=kf/Kc',i);
            for j=1:length(species)
                if (nur(j,i)>0)
                    for k=1:nur(j,i)
                        fprintf(fid, '*Y(%g)',j);
                    end
                end
            end
            fprintf(fid, '\n');
        end
        
        fprintf(fid, '\n');
        continue;
    end
    %%%%%troe Pressuredependence%%%%%start

    fprintf(fid, 'C     lindemann approach\n');
    abeLow=strread(Low{i},'%s','delimiter',' ');
    fprintf(fid, '      k0=(T**(%s))*EXP(%e - (%e)/T)\n',abeLow{2},log(str2num(abeLow{1})),str2num(abeLow{3})*41840000/RU);
    fprintf(fid, '      kinf=KFARRAY(%g)\n',i);


    fprintf(fid, '      Z=sumY\n');
    for j=1:2:(length(Meff{i})-1)
        I_species=find(strcmp(species,Meff{i}{j}));
        fprintf(fid, '      Z=Z+Y(%g)*(%g)\n',I_species,str2num(Meff{i}{j+1})-1);
    end
    fprintf(fid, '      Pr=k0*Z/kinf\n');

    %%%if after delimited troeform legnth =3 then do this if 4 then do this
    if (~isempty(Troe{i}))
        fprintf(fid, 'C     Troe form\n');
        aTTT=strread(Troe{i},'%s','delimiter',' '); 
        fprintf(fid, '      a=%s\n',aTTT{1});%%alpha
        fprintf(fid, '      T3=%s\n',aTTT{2});%%T***
        fprintf(fid, '      T1=%s\n',aTTT{3});%%T*
        switch length(aTTT)
            case 3
            	fprintf(fid, '      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1);\n');%%T** not included
            case 4
            	fprintf(fid, '      T2=%s\n',aTTT{4});%%T**
            	fprintf(fid, '      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);\n');
            otherwise
                disp('not valid troe input');
                return;
        end
        fprintf(fid, '      c = -0.4 - 0.67*LOG10(Fcent)\n');
        fprintf(fid, '      n = 0.75 - 1.27*LOG10(Fcent)\n');
        
        fprintf(fid, '      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)\n     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))\n');
    else
        fprintf(fid, '      F=1\n');
    end
    fprintf(fid, '      kf=kinf*(Pr/(1+Pr))*F\n');
    fprintf(fid, '      FWDK(%g)=kf',i);
    for j=1:length(species)
        if (nuf(j,i)>0)
            for k=1:nuf(j,i)
                fprintf(fid, '*Y(%g)',j);
            end
        end
    end
    fprintf(fid, '\n');

    if (strcmp(Rtype{i},'norm'))
        fprintf(fid, '      X=0\n');
        for j=1:length(species)
            if (nu(j,i)~=0);
                fprintf(fid, '      X=X+g(%g)*(%g)\n',j,nu(j,i));
            end
        end
        fprintf(fid, '      Kc=EXP(-X/T)');
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
        fprintf(fid, '\n');
        
        fprintf(fid, '      REVK(%g)=kf/Kc',i);
        for j=1:length(species)
            if (nur(j,i)>0)
                for k=1:nur(j,i)
                    fprintf(fid, '*Y(%g)',j);
                end
            end
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, '\n');
    
    %%%%%troe Pressuredependence%%%%%end
   

end
fprintf(fid, '      RETURN\n      END\n');
fclose(fid);

end

%things to remember:
%normal or forward reacitons, this info is inputed
%rev from princetonmech, have to add that keyword in
%pressure depenedence of species w/o troe, this info is inputed

%if troe or pressure dependence kf is kinf

%now need to look at how reg/Pdep/troe differs and how can setup this
    %code
