function [PRC1,PRC2,prc_T1,prc_T2]= getprc(t)
dt = 1e-6;
step = 15;
I = zeros(size(t));
V = hh(t,I,20);
refindex=find(findspike(V));
prc_N = refindex(end-1)-refindex(end-2);
prc_T1 = prc_N*dt;
PRC1 = zeros(1,prc_N);
for i=1:step:prc_N
    I = zeros(size(t));
    I(i+prc_N) = 2500;
    V = hh(t,I,20);
    shiftindex=find(findspike(V));
    PRC1(i:(i+step-1))= (shiftindex(end)-refindex(end))/prc_N*2*pi;
end
PRC1(PRC1>pi) = PRC1(PRC1>pi) - 2*pi;
PRC1(PRC1<-pi) = PRC1(PRC1<-pi) + 2*pi;
%plot((1:prc_N)*dt,PRC1);


V = hh(t,I,10);
refindex=find(findspike(V));
prc_N = refindex(end-1)-refindex(end-2);
prc_T2 = prc_N*dt;
PRC2 = zeros(1,prc_N);
for i=1:step:prc_N
    I = zeros(size(t));
    I(i+prc_N) = 2500;
    V = hh(t,I,10);
    shiftindex=find(findspike(V));
    PRC2(i:(i+step-1))= (shiftindex(end)-refindex(end))/prc_N*2*pi;
end
PRC2(PRC2>pi) = PRC2(PRC2>pi) - 2*pi;
PRC2(PRC2<-pi) = PRC2(PRC2<-pi) + 2*pi;
%plot((1:prc_N)*dt,PRC2);
end

