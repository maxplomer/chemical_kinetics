clc;clear;
T=linspace(300,5000,100);

for i=1:length(T)
    if T(i)<1000
        Cp_over_R(i) = (2.54205966E+00)  +  (-2.75506191E-05)*T(i)  +  (-3.10280335E-09)*T(i)^2  +  (4.55106742E-12)*T(i)^3  +  (-4.36805150E-16)*T(i)^4;
    else
        Cp_over_R(i) = (2.94642878E+00)  +  (-1.63816649E-03)*T(i)  +  (2.42103170E-06)*T(i)^2  +  (-1.60284319E-09)*T(i)^3  +  (3.89069636E-13)*T(i)^4;
    end
end

plot(T,Cp_over_R)
ylabel('Cp/R');
xlabel('T(K)')
xlim([300 5000]);

%might want to just calculate all cp/R's then linear interpolate from them