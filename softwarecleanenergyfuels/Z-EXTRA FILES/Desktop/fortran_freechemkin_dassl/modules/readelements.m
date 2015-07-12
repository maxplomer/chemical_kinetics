function readelements

global fid elements;
elements=[];
while 1
    tline = fgetl(fid);
    tline=textscan(tline,'%s');
    tline=tline{1};
    if (isempty(tline))             %deal with blank lines
       continue; 
    end
    if(strcmp(tline{1}(1),'!'))     %deal with commented lines
        continue;
    end
    if (strcmp(tline{1},'END'))     %deal with end keyword
        break;
    end
    elements=cat(1, elements, tline);
end
elements(1)=[];
end
