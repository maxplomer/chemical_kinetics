
x=-1:0.1:1;
y=-1:0.1:1;
for i=1:length(x)
    for j=1:length(y)
        f(j,i)=x(i)+y(j);
    end
end
ang=0:0.01:2*pi;
gx=cos(ang);
gy=sin(ang);
mesh(x,y,f)
alpha(.4)
hold on;
plot(gx,gy)
xlabel('x')
ylabel('y')

