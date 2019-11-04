t1=0:0.1:100;
t2=0:0.01:100;
t=(0:0.01:200)';
QA=zeros(20001,1);
FlagA = ones(680,1);
FlagB = zeros(1000,1);
FlagC = [FlagA ; FlagB];
for i=1:20001
    j=mod(i,1680);
    if j==0
        j=1680;
    end
    FlagA(i)=FlagC(j);
end
load QB;
QB=interp1(t1,QB,t2,'linear')';  
QB=[ zeros(4500,1);QB]
QB = [QB(1:10000) ; QB(1:10001)];
PA = zeros(20001,1);
PB = zeros(20001,1);
RhoA = zeros(20001,1);
RhoB = zeros(20001,1);

C=0.85;
A=pi*0.7*0.7;
V=pi*5*5*500;
PA(1) = 160;
PB(1) = 100;
RhoA(1) = fun_PreToDen(PA(1));
RhoB(1) = fun_PreToDen(PB(1));
QA(1) = C * A * sqrt( 2*(PA(1)-PB(1))/RhoA(1) ) ;

for i=2:20001
    PA(i) = PA(i-1);
    RhoA(i) = fun_PreToDen(PA(i));  
    RhoB(i) = RhoB(i-1) + (tt(i) - tt(i-1)) * ( QA(i-1) * RhoA(i-1) - QB(i-1) * RhoB(i-1) ) / V;   
    PB(i) = fun_DenToPre(RhoB(i));  
    if FlagA(i)==1
        QA(i) = C * A * sqrt( 2*(PA(i)-PB(i))/RhoA(i) ) ;
    else
        QA(i) = 0;
    end 
end




