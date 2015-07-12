function plotoutput
clc;clear;

[tS,y1S,y2S,y3S]=textread('output.txt','%s %s %s %s');

for i=1:length(tS)
	t(i)=str2num(tS{i});
	y1(i)=str2num(y1S{i});
	y2(i)=str2num(y2S{i});
	y3(i)=str2num(y3S{i});
end
% plot the results

figure

plot(t, y1)
%hold on;
%plot(t, y2,'.')
%hold on;
%plot(t, y3,'x')






end
