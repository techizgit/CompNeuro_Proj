clc;clear;
Omega = 2*pi*20;
dt = 1e-6;
M = 5;
T = 2*pi*M/Omega;
t = (-T/2):dt :(T/2);
t2 = 0:dt:T;
H = [zeros(1,125000),ones(1,125001)];
Hd = [zeros(1,1000),ones(1,249001)];
Hd2 = [zeros(1,500),ones(1,249501)];
h1 = 10*exp(-100*t).*((100*t).^3/6-(100*t).^5/120).*H;
h2 = 10*exp(-100*t).*((100*t).^5/120).*H;
h11 = -5*exp(-200*(t2-0.001)).*(((200*(t2-0.001)).^3)/6).*Hd;
h12 = 10*exp(-200*(2*t2-0.001)).*(((200*(2*t2-0.001)).^5)/120-((200*(2*t2-0.001)).^7)/5040).*Hd2;
h21 = -h11;
h22 = -20*exp(-200*(t2-0.001)).*(((200*(t2-0.001)).^5)/120-((200*(t2-0.001)).^7)/5040).*Hd;

ut = zeros(size(t));
m = -M:M;
for i = 1:length(m)
    ut = ut + real(exp(1i*m(i)*Omega*t/M));
end
c = 10;
ut = c*ut/max(ut);

uh1 = conv(ut,h1)*dt;uh1=uh1(125000:375000);
uh2 = conv(ut,h2)*dt;uh2=uh2(125000:375000);
%% 

BSG1 = zeros(1,length(t)); 
BSG1(1) = 50;
BSG2 = zeros(1,length(t)); 
BSG2(1) = 50;
v1 = zeros(1,length(t)); 
v1(1) = uh1(1);
v2 = zeros(1,length(t)); 
v2(1) = uh2(1);
x1 = [0 0 1];
x2 = [0 0 1];
g = [120 36 0.3];
spikenumber1=0;
spikenumber2=0;

dt = 1000*dt;
for i=2:length(t)
    v1(i) = uh1(i);
    v2(i) = uh2(i);
    if spikenumber1>0
        for j=1:spikenumber1
            v1(i)=v1(i)+h11(i-spikeindex1(j));
            v2(i)=v2(i)+h12(i-spikeindex1(j));
        end
    end
    if spikenumber2>0
        for j=1:spikenumber2
            v1(i)=v1(i)+h21(i-spikeindex2(j));
            v2(i)=v2(i)+h22(i-spikeindex2(j));
        end
    end
    
    a = exp(-(BSG1(i-1)+55)/10)-1;
    if a == 0
        dn = (1-x1(1)) * 0.1 - x1(1) * (0.125*exp(-(BSG1(i-1)+65)/80));
    else
        dn = (1-x1(1)) * (-0.01*(BSG1(i-1)+55)/a) - x1(1) * (0.125*exp(-(BSG1(i-1)+65)/80));
    end
    a = exp(-(BSG1(i-1)+40)/10)-1;
    if a == 0
        dm = (1-x1(2)) - x1(2) * (4*exp(-(BSG1(i-1)+65)/18));
    else
        dm = (1-x1(2)) * (-0.1*(BSG1(i-1)+40)/a) - x1(2) * (4*exp(-(BSG1(i-1)+65)/18));
    end
    dh = (1-x1(3)) * (0.07*exp(-(BSG1(i-1)+65)/20)) - x1(3) / (exp(-(BSG1(i-1)+35)/10)+1);
    x1 = x1 + dt*[dn dm dh];
    BSG1(i) = BSG1(i-1)+dt*(- g(1)*x1(2).^3*x1(3)*(BSG1(i-1)-50) - g(2)*x1(1).^4*(BSG1(i-1)+77) - g(3)*(BSG1(i-1)+54.387)+20+v1(i));
    if i>2
        if (BSG1(i)<BSG1(i-1))&&(BSG1(i-2)<BSG1(i-1))
            spikenumber1=spikenumber1+1;
            spikeindex1(spikenumber1)=i-1;
        end
    end
    
    a = exp(-(BSG2(i-1)+55)/10)-1;
    if a == 0
        dn = (1-x2(1)) * 0.1 - x2(1) * (0.125*exp(-(BSG2(i-1)+65)/80));
    else
        dn = (1-x2(1)) * (-0.01*(BSG2(i-1)+55)/a) - x2(1) * (0.125*exp(-(BSG2(i-1)+65)/80));
    end
    a = exp(-(BSG2(i-1)+40)/10)-1;
    if a == 0
        dm = (1-x2(2)) - x(2) * (4*exp(-(BSG2(i-1)+65)/18));
    else
        dm = (1-x2(2)) * (-0.1*(BSG2(i-1)+40)/a) - x2(2) * (4*exp(-(BSG2(i-1)+65)/18));
    end
    dh = (1-x2(3)) * (0.07*exp(-(BSG2(i-1)+65)/20)) - x2(3) / (exp(-(BSG2(i-1)+35)/10)+1);
    x2 = x2 + dt*[dn dm dh];
    BSG2(i) = BSG2(i-1)+dt*(- g(1)*x2(2).^3*x2(3)*(BSG2(i-1)-50) - g(2)*x2(1).^4*(BSG2(i-1)+77) - g(3)*(BSG2(i-1)+54.387)+10+v2(i));
    if i>2
        if (BSG2(i)<BSG2(i-1))&&(BSG2(i-2)<BSG2(i-1))
            spikenumber2=spikenumber2+1;
            spikeindex2(spikenumber2)=i-1;
        end
    end
end
dt = 1e-6;

figure();
subplot(211);plot(t,BSG1);
hold on;
plot(t(spikeindex1),BSG1(spikeindex1), '*');ylabel('HHN1 Output');
hold off;xlabel('t(s)');
subplot(212);plot(t,BSG2);
hold on;
plot(t(spikeindex2),BSG2(spikeindex2), '*');ylabel('HHN2 Output');
hold off;xlabel('t(s)');
%% 


%[PRC1,PRC2,prc_T1,prc_T2]=getprc(t);


I = zeros(size(t));
V = hh(t,I,20);
refindex=find(findspike(V));
prc_N = refindex(end-1)-refindex(end-2);
prc_T1 = prc_N*dt;
V = hh(t,I,10);
refindex=find(findspike(V));
prc_N = refindex(end-1)-refindex(end-2);
prc_T2 = prc_N*dt;
load('PRC1.mat','PRC1');PRC1=PRC1*2;
load('PRC2.mat','PRC2');PRC2=PRC2*2;

figure();
subplot(211);plot((1:length(PRC1))*dt,PRC1);ylabel('HHN1 PRC');
xlabel('t(s)');
subplot(212);plot((1:length(PRC2))*dt,PRC2);ylabel('HHN2 PRC');
xlabel('t(s)');
%% 
v1 = zeros(1,length(t)); 
v1(1) = uh1(1);
v2 = zeros(1,length(t)); 
v2(1) = uh2(1);
intv1=0;
intv2=0;
PIF1 = zeros(size(t));
PIF2 = zeros(size(t));
i_prc1 =1;
i_prc2 =1;
pspikenumber1=0;
pspikenumber2=0;
PRCrep1 = [PRC1 PRC1];
PRCrep2 = [PRC2 PRC2];
th1 = prc_T1*20;
th2 = prc_T2*10;
for i=1:length(t)
    v1(i) = uh1(i);
    v2(i) = uh2(i);
    if pspikenumber1>0
        for j=1:pspikenumber1
            v1(i)=v1(i)+h11(i-pspikeindex1(j));
            v2(i)=v2(i)+h12(i-pspikeindex1(j));
        end
    end
    if pspikenumber2>0
        for j=1:pspikenumber2
            v1(i)=v1(i)+h21(i-pspikeindex2(j));
            v2(i)=v2(i)+h22(i-pspikeindex2(j));
        end
    end
    
    intv1=intv1+dt*(20+v1(i)*PRCrep1(i_prc1));
    if intv1>th1
        intv1 = intv1 - th1;
        i_prc1 = 1;
        pspikenumber1 = pspikenumber1+1;
        pspikeindex1(pspikenumber1)=i;
    else
        i_prc1 = i_prc1+1;
    end 
    PIF1(i)=intv1;
    
    intv2=intv2+dt*(10+v2(i)*PRCrep2(i_prc2));
    if intv2>th2
        intv2 = intv2 - th2;
        i_prc2 = 1;
        pspikenumber2 = pspikenumber2+1;
        pspikeindex2(pspikenumber2)=i;
    else
        i_prc2 = i_prc2+1;
    end
    PIF2(i)=intv2;
end

figure();
subplot(211);plot(t,PIF1);
hold on;
plot(t(pspikeindex1),PIF1(pspikeindex1), '*');ylabel('PIF1 Output');
hold off;xlabel('t(s)');
subplot(212);plot(t,PIF2);
hold on;
plot(t(pspikeindex2),PIF2(pspikeindex2), '*');ylabel('PIF2 Output');
hold off;xlabel('t(s)');
%% 
spiketime11=spikeindex1(2:end)*dt;
spiketime22=spikeindex2(2:end)*dt;
pspiketime11=pspikeindex1*dt;
pspiketime22=pspikeindex2(1:end-1)*dt;
figure();
subplot(221);plot((1:21),spiketime11,'--*',(1:21),pspiketime11,'--*');
ylabel('HHN1 Spike Time & PIF1 Spike Time');
subplot(222);plot((1:16),spiketime22,'--*',(1:16),pspiketime22,'--*');
ylabel('HHN2 Spike Time & PIF2 Spike Time');
subplot(223);plot((1:21),spiketime11-pspiketime11,'--*');
ylabel('HHN1-PIF1 Spike Time Difference');
subplot(224);plot((1:16),spiketime22-pspiketime22,'--*');
ylabel('HHN2-PIF2 Time Difference');


%% 

[q1,q2]=qcal(t,pspikeindex1, pspikeindex2, PRCrep1, PRCrep2, th1, th2, 20, 10);
%% 

gt=gfunc(t);
g11=conv(conv(h1,fliplr(h1))*dt,gt*dt);
g12=conv(conv(h1,fliplr(h2))*dt,gt*dt);
g21=conv(conv(h2,fliplr(h1))*dt,gt*dt);
g22=conv(conv(h2,fliplr(h2))*dt,gt*dt);
G11 = zeros(pspikenumber1-1,pspikenumber1-1);
for k=1:pspikenumber1-1
    for l=1:pspikenumber1-1
        G11(k,l)=0;
        for s=pspikeindex1(k):pspikeindex1(k+1)
            G11(k,l)=G11(k,l)+g11(s-ceil((pspikeindex1(l+1)+pspikeindex1(l))/2)+375000)*PRCrep1(s-pspikeindex1(k)+1)*dt;
        end
    end
end
G12 = zeros(pspikenumber1-1,pspikenumber2-1);
for k=1:pspikenumber1-1
    for l=1:pspikenumber2-1
        G12(k,l)=0;
        for s=pspikeindex1(k):pspikeindex1(k+1)
            G12(k,l)=G12(k,l)+g12(s-ceil((pspikeindex2(l+1)+pspikeindex2(l))/2)+375000)*PRCrep1(s-pspikeindex1(k)+1)*dt;
        end
    end
end
G21 = zeros(pspikenumber2-1,pspikenumber1-1);
for k=1:pspikenumber2-1
    for l=1:pspikenumber1-1
        G21(k,l)=0;
        for s=pspikeindex2(k):pspikeindex2(k+1)
            G21(k,l)=G21(k,l)+g21(s-ceil((pspikeindex1(l+1)+pspikeindex1(l))/2)+375000)*PRCrep2(s-pspikeindex2(k)+1)*dt;
        end
    end
end
G22 = zeros(pspikenumber2-1,pspikenumber2-1);
for k=1:pspikenumber2-1
    for l=1:pspikenumber2-1
        G22(k,l)=0;
        for s=pspikeindex2(k):pspikeindex2(k+1)
            G22(k,l)=G22(k,l)+g22(s-ceil((pspikeindex2(l+1)+pspikeindex2(l))/2)+375000)*PRCrep2(s-pspikeindex2(k)+1)*dt;
        end
    end
end
%%
%test1=zeros(1,pspikenumber1-1);
%test2=zeros(1,pspikenumber2-1);
%%%%%%%%%%%%%%%%%%%%
% 
% for i=1:pspikenumber1-1
%     intv1=0;
%     ttt=1;
%     for s = pspikeindex1(i):pspikeindex1(i+1)
%         intv1=intv1+uh1(s)*PRCrep1(ttt)*dt;
%         ttt=ttt+1;
%     end
%     test1(i)=intv1;
% end
% 
% for i=1:pspikenumber2-1
%     intv2=0;
%     ttt=1;
%     for s = pspikeindex2(i):pspikeindex2(i+1)
%         intv2=intv2+uh2(s)*PRCrep2(ttt)*dt;
%         ttt=ttt+1;
%     end
%     test2(i)=intv2;
% end
%%%%%%%%%%%%%%%%%%%%%
%%
G = [G11,G12;G21,G22];
q = [q1,q2]';
c = pinv(G,1e-8)*q;
%tq = [test1,test2]';
%tc = pinv(G,1e-8)*tq;

c1 = c(1:pspikenumber1-1);
c2 = c(pspikenumber1:end);
%%
phi1=conv(fliplr(h1),gt*dt);
phi2=conv(fliplr(h2),gt*dt);
ut_rec = zeros(size(t));

for i = 1:length(t)
    for j = 1:length(c1)
        ut_rec(i)=ut_rec(i)+c1(j)*phi1(i-ceil((pspikeindex1(j+1)+pspikeindex1(j))/2)+250000);
    end
    for j = 1:length(c2)
        ut_rec(i)=ut_rec(i)+c2(j)*phi2(i-ceil((pspikeindex2(j+1)+pspikeindex2(j))/2)+250000);
    end
end
%%
figure();
plot(t,ut_rec,t,ut);xlabel('t(s)');ylabel('Input and Recovered Input (c=10)');
MSE(1) = immse(ut/max(ut),ut_rec/max(ut_rec));
MSE(2) = immse(ut,ut_rec);

error = 10*log10(abs(ut_rec-ut));
figure();plot(t,error);xlabel('t(s)');
ylabel('Error Magnitude(dB)');

