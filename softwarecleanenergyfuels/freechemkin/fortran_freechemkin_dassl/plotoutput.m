function plotoutput

[t,y1,y2,y3,y4]=textread('output.txt','%f %f %f %f %f');


% plot the results

figure(1)
plot(t, y1)
xlabel('Time (sec)');
ylabel('Temperature (K)');

figure(2)
plot(t, y2)
xlabel('Time (sec)');
ylabel('Mass Fraction CH4 (g/g)');

figure(3)
plot(t, y3)
xlabel('Time (sec)');
ylabel('Mass Fraction H radical (g/g)');

figure(4)
plot(t, y4)
xlabel('Time (sec)');
ylabel('Mass Fraction CO (g/g)');





end
