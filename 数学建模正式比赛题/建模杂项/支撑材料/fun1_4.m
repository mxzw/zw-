t1=0:0.1:100;
t2=0:0.01:100;
t=(0:0.01:500)';
QA=zeros(50001,1);
FlagA = ones(230,1);
FlagB = zeros(1000,1);
FlagC = [FlagA ; FlagB];
for i=1:50001
    j=mod(i,1230);
    if j==0
        j=1230;
    end
    GateA(i)=FlagC(j);
end
load QB;
QB1=interp1(t1,QB,t2,'linear')';  
QB2=[ zeros(4500,1);QB1]
QB2 = [QB2(1:10000) ; QB2(1:10000); QB2(1:10000); QB2(1:10000); QB2(1:10001)];
PA = zeros(50001,1);
PB = zeros(50001,1);
RhoA = zeros(50001,1);
RhoB = zeros(50001,1);
c=0.85;
a=pi*0.7*0.7;
v=pi*5*5*500;
PA(1) = 160;
PB(1) = 100;
RhoA(1) = fun_PreToDen(PA(1));
RhoB(1) = fun_PreToDen(PB(1));
QA(1) = c * a * sqrt( 2*(PA(1)-PB(1))/RhoA(1) ) ;
for i=2:50001
    PA(i) = PA(i-1);
    RhoA(i) = fun_PreToDen(PA(i));
    RhoB(i) = RhoB(i-1) + (t(i) - t(i-1)) * ( QA(i-1) * RhoA(i-1) - QB2(i-1) * RhoB(i-1) ) / v;
    PB(i) = fun_DenToPre(RhoB(i)); 
    if GateA(i)==1
        QA(i) = c * a * sqrt( 2*(PA(i)-PB(i))/RhoA(i) ) ;
    else
        QA(i) = 0;
    end 
end




