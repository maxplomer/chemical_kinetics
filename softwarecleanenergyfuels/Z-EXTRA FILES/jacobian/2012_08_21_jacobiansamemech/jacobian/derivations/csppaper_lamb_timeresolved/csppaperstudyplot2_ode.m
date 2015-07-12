function csppaperstudyplot2_ode
clc;clear;
format longe;
R=3;
KK=3;
k1=1;k2=2;k3=3;

y0(1)=1;
y0(2)=0;
y0(3)=0;

S1=[-1 1 0];
S2=[0 -1 1];
S3=[0 1 -1];
S=[S1 ; S2 ; S3];



J=[ dot(S1,DF1DY)   dot(S1,DF2DY)    dot(S1,DF3DY) ;
    dot(S2,DF1DY)   dot(S2,DF2DY)    dot(S2,DF3DY) ;
    dot(S3,DF1DY)   dot(S3,DF2DY)    dot(S3,DF3DY) ];


[B,~] = eig(J);
mu=eig(J);
A=B^-1;

% test equation 3.5
total=0;
for m=1:R
    for r=1:R
        total1 =mu(m)*B(r,m);
        total2=0;
        for k=1:R
            total2=total2 +B(k,m)*J(r,k);
        end
        total1;
        total2;
        total=total-total1+total2;
    end
end
%total%=0




%time span 0 to 10 seconds
tspan = [0 10]; %

% setup and run the matlab integration solver
options = odeset('RelTol', 1.0E-6, 'AbsTol', 1.0E-9);
[t ,y] = ode15s(@ignfun, tspan, y0, options);


a=zeros(R,KK);
for m=1:R
    for r=1:R
        a(m,:)=a(m,:) + S(r,:)*A(m,r); 
    end
end




for i=1:length(t)
    f=zeros(1,R);
    F=getF(y(i,:));
    for m=1:R
        for r=1:R
            f(m)=f(m) + B(r,m)*F(r);
        end
    end
    K=zeros(R,KK);
    for m=1:R
        for s=1:KK
            term1=0;
            for k=1:R
               term1=term1+abs(a(k,s)*f(k)); 
            end
            K(m,s)=abs(a(m,s)*f(m))/term1;
        end
    end

    
    
    K_R1_S1(i)=K(1,1);
    K_R1_S2(i)=K(1,2);
    K_R1_S3(i)=K(1,3);
    
    
    K_R2_S1(i)=K(2,1);
    K_R2_S2(i)=K(2,2);
    K_R2_S3(i)=K(2,3);
    
    K_R3_S1(i)=K(3,1);
    K_R3_S2(i)=K(3,2);
    K_R3_S3(i)=K(3,3);
      
    
end

figure
plot(t,y(:,1),t,y(:,2),'x',t,y(:,3))



figure
plot(t,K_R1_S1)
figure
plot(t,K_R1_S2)
figure
plot(t,K_R1_S3)

figure
plot(t,K_R2_S1)
figure
plot(t,K_R2_S2)
figure
plot(t,K_R2_S3)

figure
plot(t,K_R3_S1)
figure
plot(t,K_R3_S2)
figure
plot(t,K_R3_S3)

    function dydt = ignfun(t, y)
        
        %dydt = S1*F1(y)+S2*F2(y);
        
%         dydt = 0;
%         F=getF(y);
%         for r=1:R
%             dydt = dydt + S(r,:)*F(r);
%         end
%         dydt=dydt';
        


        a=zeros(R,KK);
        for m=1:R
            for r=1:R
                a(m,:)=a(m,:) + S(r,:)*A(m,r); 
            end
        end
        f=zeros(1,R);
        F=getF(y);
        for m=1:R
            for r=1:R
                f(m)=f(m) + B(r,m)*F(r);
            end
        end
        dydt = 0;
        for m=1:R
            dydt = dydt + a(m,:)*f(m);
        end
        dydt=dydt';


%         F=getF(y);
%         for k=1:KK
%             dydt(k) = sum(S(k,:)'.*F');
%         end
%         dydt=dydt';
    end




    function [F]=getF(y)
        F=[F1(y) F2(y) F3(y)];
    end
    function [F1]=F1(y)
        F1=k1*y(1);
    end
    function [F2]=F2(y)
        F2=k2*y(2);
    end
    function [F3]=F3(y)
        F3=k3*y(3);
    end



    function [DFDY]=getDFDY
        DFDY=[DF1DY ; DF2DY  ; DF3DY];
    end
    function [DF1DY]=DF1DY
        DF1DY=[k1 0 0];
    end
    function [DF2DY]=DF2DY
        DF2DY=[0 k2 0];
    end
    function [DF3DY]=DF3DY
        DF3DY=[0 0 k3];
    end




end
