function  reactor_3dplotcode_plot


load data.mat

figure(1)
mesh(f,R,ignitiondelay)

ylabel('R, compression ratio')
xlabel('f, frequency (Hz)')
zlabel('time of Tmax (sec) (Note: 0 is no ignition)')

figure(2)
mesh(f,R,ignitiondelay2)

ylabel('R, compression ratio')
xlabel('f, frequency (Hz)')
zlabel('Dt from half fuel mass fraction to Tmax (sec) (Note: 0 is no ignition)')

figure(3)
mesh(f,R,ignition)

ylabel('R, compression ratio')
xlabel('f, frequency')
zlabel('0=no ignition, 1=ignition')

figure(4)
mesh(f,R,Tmax)

ylabel('R, compression ratio')
xlabel('f, frequency')
zlabel('Max Temp (K)')

end