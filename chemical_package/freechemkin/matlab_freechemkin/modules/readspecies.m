function readspecies

global fid species;
species=[];
while 1
    tline = fgetl(fid);
    tline=strread(tline,'%s','delimiter',' ');
    if (isempty(tline))             %deal with blank lines
       continue; 
    end
    if(strcmp(tline{1}(1),'!'))     %deal with commented lines
        continue;
    end
    if (strcmp(tline{1},'END'))     %deal with end keyword
        break;
    end
    species=cat(1, species, tline);
end

end