t1=0:0.1:100;
t2=0:0.01:100;
t=(0:0.01:1000)';
QA=zeros(100001,1);
FlagA = ones(130,1);
FlagB = zeros(1000,1);

FlagC = [FlagA ; FlagB];
for I=1:100001
    J=mod(I,1130);
    if J==0
        J=1130;
    end
    FlagA(I)=FlagC(J);
end
load QB
QB1=interp1(t1,QB,t2,'linear')';  
QB2=[ zeros(4500,1);QB1]
QB2 = [QB2(1:10000) ; QB2(1:10000); QB2(1:10000); QB2(1:10000); QB2(1:10000) ; QB2(1:10000); QB2(1:10000); QB2(1:10000); QB2(1:10000); QB2(1:10001)];
PA = zeros(100001,1);
PB = zeros(100001,1);
RhoA = zeros(100001,1);
RhoB = zeros(100001,1);

C=0.85;
A=pi*0.7*0.7;
V=pi*5*5*500;

PA(1) = 160;
PB(1) = 100;

RhoA(1) = fun_PreToDen(PA(1));
RhoB(1) = fun_PreToDen(PB(1));

QA(1) = C * A * sqrt( 2*(PA(1)-PB(1))/RhoA(1) ) ;

%%
for I=2:100001
    PA(I) = PA(I-1);
    RhoA(I) = fun_PreToDen(PA(I));
    
    RhoB(I) = RhoB(I-1) + (t(I) - t(I-1)) * ( QA(I-1) * RhoA(I-1) - QB2(I-1) * RhoB(I-1) ) / V;
    
    PB(I) = fun_DenToPre(RhoB(I));
    
    if FlagA(I)==1
        QA(I) = C * A * sqrt( 2*(PA(I)-PB(I))/RhoA(I) ) ;
    else
        QA(I) = 0;
    end
    
end





