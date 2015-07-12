function [DfwdkDC,DrevkDC]=getDkfkrDC21(T,C,j)
%reaction number
I=21;


RU=83145100; %erg/(mol*K)
[A,B,E]=getabe; %-units of E are cal/mol , 1 calorie = 41 840 000 erg
                %-units of A are cgs (cm, sec, K, mole), the exact units 
                % depending on the reaction.
[nuf,nur]=getnu;                            
kf=A(I)*T^B(I)*exp(-E(I)*41840000/(RU*T));

%C is molar concentration c=X*P/(RU*T)
DfwdkDC=kf;
for i=1:9
   if i==j
       DfwdkDC=DfwdkDC*nuf(i,I)*C(i)^(nuf(i,I)-1);
   else
       DfwdkDC=DfwdkDC*C(i)^nuf(i,I);
   end
end


%Kp(T)=exp{-[sum i=1:N of (vi''-vi')*mu0,i(T)]/(RU*T)}

g=getg(T);
Kp=exp(   -dot(nur(:,I)-nuf(:,I),g)  /  (RU*T)  );

%kf/kb=Kc=Kp/(RU*T)^[sum i=1:N of vi'' - vi']
PATM=1.01325D6;
Kc=Kp/(RU*T/PATM)^sum(nur(:,I)-nuf(:,I));
kb=kf/Kc;

%C is molar concentration c=X*P/(RU*T)
DrevkDC=kb;
for i=1:9
    if i==j
        DrevkDC=DrevkDC*nur(i,I)*C(i)^(nur(i,I)-1);
    else
        DrevkDC=DrevkDC*C(i)^nur(i,I);
    end
end

end