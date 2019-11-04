t1=(0:0.1:200)';
QB=zeros(2001,1);
flagA = ones(1001,1);
for i=29:1001
    flagA(i)=0;
end

flagA = [flagA(1:1000) ; flagA]

load QB;
QB=[ zeros(2,1);QB]
QB = [QB(1:1000) ; QB(1:1001)];
PA = zeros(2001,1);
PB = zeros(2001,1);

RhoA = zeros(2001,1);
RhoB = zeros(2001,1);
C=0.85;
A=pi*0.7*0.7;
V=pi*5*5*500;
PA(1) = 160;
PB(1) = 100;
RhoA(1) = fun_PreToDen(PA(1));
RhoB(1) = fun_PreToDen(PB(1));
QB(1) = C * A * sqrt( 2*(PA(1)-PB(1))/RhoA(1) ) ;
for i=2:2001
    PA(i) = PA(i-1);
    RhoA(i) = fun_PreToDen(PA(i));
    RhoB(i) = RhoB(i-1) + (t1(i) - t1(i-1)) * ( QB(i-1) * RhoA(i-1) - QB(i-1) * RhoB(i-1) ) / V;
    PB(i) = fun_DenToPre(RhoB(i));
    if flagA(i)==1
        QB(i) = C * A * sqrt( 2*(PA(i)-PB(i))/RhoA(i) ) ;
    else
        QB(i) = 0;
    end
end




