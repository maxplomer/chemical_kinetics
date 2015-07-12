function plotoutput
fid = fopen('output.txt','r');
z=0;
t=[];
while 1
	z=z+1;
	
	tline = fgetl(fid);
	if isempty(tline)
		break;
	end

	t(z)=str2num(tline(10:20));
	tline=tline(28:length(tline));
	tline=textscan(tline,'%s');
	tline=tline{1};
	for i=1:length(tline)
		y(z,i)=str2num(tline{i});
	end
	
end





% plot the results
 figure

plot(t, y(:, 1))
hold on;
plot(t, y(:, 2),'.')
hold on;
plot(t, y(:, 3),'x')






end
