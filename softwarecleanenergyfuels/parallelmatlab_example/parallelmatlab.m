clc;clear;
tic
for i=1:4
    wastetime;
end
toc
tic
matlabpool(4)
toc
tic
parfor i=1:4
    wastetime;
end
toc
tic
matlabpool('close');
toc