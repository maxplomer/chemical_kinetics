function cp=getDcpDT(T)
%cp=cv + RU
%units: erg/(mol*K)     erg=g*cm^2/s^2=1e-7 J
RU=83145100; %erg/(mol*K)


% k_names = 
% 
%     'H2'
%     'O2'
%     'H2O'
%     'H'
%     'O'
%     'OH'
%     'HO2'
%     'H2O2'
%     'N2'



%     'H2'
A(:,1)=[ 2.99142337E+00 7.00064411E-04 -5.63382869E-08 -9.23157818E-12 1.58275179E-15    ...
-8.35033997E+02 -1.35511017E+00 3.29812431E+00 8.24944174E-04 -8.14301529E-07    ...
-9.47543433E-11 4.13487224E-13 -1.01252087E+03 -3.29409409E+00];                
%     'O2'
A(:,2)=[ 3.69757819E+00 6.13519689E-04 -1.25884199E-07 1.77528148E-11 -1.13643531E-15    ...
-1.23393018E+03 3.18916559E+00 3.21293640E+00 1.12748635E-03 -5.75615047E-07    ...
 1.31387723E-09 -8.76855392E-13 -1.00524902E+03 6.03473759E+00];
%     'H2O'
A(:,3)=[ 2.67214561E+00 3.05629289E-03 -8.73026011E-07 1.20099639E-10 -6.39161787E-15    ...
-2.98992090E+04 6.86281681E+00 3.38684249E+00 3.47498246E-03 -6.35469633E-06    ...
 6.96858127E-09 -2.50658847E-12 -3.02081133E+04 2.59023285E+00];
%     'H'
A(:,4)=[2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    ...
 2.54716270E+04 -4.60117638E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    ...
 0.00000000E+00 0.00000000E+00 2.54716270E+04 -4.60117608E-01];
%     'O'
A(:,5)=[ 2.54205966E+00 -2.75506191E-05 -3.10280335E-09 4.55106742E-12 -4.36805150E-16    ...
 2.92308027E+04 4.92030811E+00 2.94642878E+00 -1.63816649E-03 2.42103170E-06    ...
-1.60284319E-09 3.89069636E-13 2.91476445E+04 2.96399498E+00];
%     'OH'
A(:,6)=[ 2.86472886E+00 1.05650448E-03 -2.59082758E-07 3.05218674E-11 -1.33195876E-15    ...
 3.68362875E+03 5.70164073E+00 4.12530561E+00 -3.22544939E-03 6.52764691E-06    ...
-5.79853643E-09 2.06237379E-12 3.34630913E+03 -6.90432960E-01 ];
%     'HO2'
A(:,7)=[ 4.01721090E+00 2.23982013E-03 -6.33658150E-07 1.14246370E-10 -1.07908535E-14    ...
 1.11856713E+02 3.78510215E+00 4.30179801E+00 -4.74912051E-03 2.11582891E-05    ...
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 ];
%     'H2O2'
A(:,8)=[4.57316685E+00 4.33613639E-03 -1.47468882E-06 2.34890357E-10 -1.43165356E-14    ...
-1.80069609E+04 5.01136959E-01 3.38875365E+00 6.56922581E-03 -1.48501258E-07    ...
-4.62580552E-09 2.47151475E-12 -1.76631465E+04 6.78536320E+00];
%     'N2'
A(:,9)=[ 0.02926640E+02 0.01487977E-01 -0.05684761E-05 0.01009704E-08 -0.06753351E-13    ...
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.01408240E-01 -0.03963222E-04    ...
 0.05641515E-07 -0.02444855E-10 -0.01020900E+05 0.03950372E+02];




if T>1000
    for i=1:9
        cp_over_R(i)=A(2,i) +2*A(3,i)*T +3*A(4,i)*T^2 +4*A(5,i)*T^3;
    end
    
else
    for i=1:9
        cp_over_R(i)=A(9,i) +2*A(10,i)*T +3*A(11,i)*T^2 +4*A(12,i)*T^3;
    end
end
cp=cp_over_R*RU;

end