function [outputcellarray]=runodecases(P) 
%P is intial condition range

hosts = [ "54.200.141.205"; 
    "54.201.215.228"; 
    "54.201.212.205"];
N=2;%number of slaves
sockets = connect (hosts);


for i=1:length(P)
    cmds{i}=sprintf("odeexample(%g)",P(i));
endfor



%returns the outputcellarray from a list of commands
i=0;
while i<length(cmds)
    for j=1:N
        if (i+j)>length(cmds)
            break;
        endif
        cmd=sprintf("send (%s, sockets(1, :))",cmds{i+j});
        reval (cmd, sockets(j+1, :));
    endfor
    for j=1:N
        if (i+j)>length(cmds)
            break;
        endif
        outputcellarray{i+j} = recv (sockets(j+1, :));
    endfor
    i=i+N;
endwhile

%%%%%%%%

scloseall (sockets);

endfunction