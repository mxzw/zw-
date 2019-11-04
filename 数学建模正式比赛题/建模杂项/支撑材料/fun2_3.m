close all;
clear all;
clc;

%%  2-1 

tt1=0:0.1:100;
tt2=0:0.01:100;
tt=(0:0.01:1000)';
%%
Omega = 0.6;
for ii=1:100001
    
    temp=mod(tt(ii)*Omega , 6.27);
    VV(ii)=fun_ThetaToV(temp);
end
VV=VV';

MM=zeros(100001,1);


%%
QA=zeros(100001,1);
QB=zeros(10001,1);
QD=zeros(4001,1);
load QB_S;
load QB_D;

QB1 = [QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ;QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ;QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ;QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ;QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ;QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ;QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ; QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ; QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4000) ; QB_S(1:3000) ; QB_S(1:3000) ; QB_D(1:4001) ;];
%%
PA = zeros(100001,1);
PB = zeros(100001,1);

RA = zeros(100001,1);
RB = zeros(100001,1);


%% 常数
C_C=0.85;
C_A=pi*0.7*0.7;
C_V=pi*5*5*500;


%% 初始条件
PA(1) = 160;
PB(1) = 100;

RA(1) = fun_PreToDen(PA(1));
RB(1) = fun_PreToDen(PB(1));

QA(1) = C_C * C_A * sqrt( 2*(PA(1)-PB(1))/RA(1) ) ;
QB(1) = C_C * QB1(1) * sqrt( 2*(PB(1)-0.1) / RB(1) ) ;

MM(1) = RA(1)*VV(1);

%%
for ii=2:100001
    
    RB(ii) = RB(ii-1) + (tt(ii) - tt(ii-1)) * ( QA(ii-1) * RA(ii-1) - QB(ii-1) * RB(ii-1) ) / C_V;
    PB(ii) = fun_DenToPre(RB(ii));
    MM(ii)=MM(ii-1);
%     if VV(ii)>VV(ii-1)
%         if PA(ii-1)>PB(ii-1)
%             MM(ii)=MM(ii-1)-QA(ii-1)*RA(ii-1)*(tt(ii) - tt(ii-1));
%         else
%             if RA(ii-1)>0.85
%                 MM(ii)=MM(ii-1);
%             else
%                 MM(ii)=0.85*VV(ii);
%             end
%         end
%     else
%         if PA(ii-1)>PB(ii-1)
%             MM(ii)=MM(ii-1)-QA(ii-1)*RA(ii-1)*(tt(ii) - tt(ii-1));
%         else
%             MM(ii)=MM(ii-1);
%         end
%     end
%     
    RA(ii) = MM(ii)/VV(ii);
    if RA(ii)<=0.85
        RA(ii)=0.85;
    end
    PA(ii) = fun_DenToPre(RA(ii));
            
    if PA(ii)>PB(ii)
        QA(ii) = C_C * C_A * sqrt( 2*(PA(ii)-PB(ii))/RA(ii) ) ;
    else
        QA(ii) = 0;
    end
    
    QB(ii) = C_C * QB1(ii) * sqrt( 2*(PB(ii)-0.1) / RB(ii) ) ;
    
end

figure
plot(tt,VV,tt,MM)
figure
plot(tt,QA)
figure
plot(tt,QB)
figure
plot(tt,PB)




