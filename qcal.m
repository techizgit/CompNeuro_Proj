%function [ u_rec, G, q ] = decode_PIF(t, spk_index1, spk_index2, omega, prc1, prc2, threshold1, threshold2, bias1, bias2)
function [q1, q2] = qcal(t,PIF_tk1, PIF_tk2, prc1, prc2, thr1, thr2, bias1, bias2)

h11 = @(t, tk) -5*exp(-200*((t-tk)-1e-3)).*(200*((t-tk)-1e-3)).^3/factorial(3).*((t-tk)>=1e-3);
h12 = @(t, tk) 10*exp(-200*(2*(t-tk)-1e-3)).*((200*(2*(t-tk)-1e-3)).^5/factorial(5) - (200*(2*(t-tk)-1e-3)).^7/factorial(7)).*(2*(t-tk)>=1e-3);
h21 = @(t, tk) 5*exp(-200*((t-tk)-1e-3)).*(200*((t-tk)-1e-3)).^3/factorial(3).*((t-tk)>=1e-3);
h22 = @(t, tk) -20*exp(-200*((t-tk)-1e-3)).*((200*((t-tk)-1e-3)).^5/factorial(5) - (200*((t-tk)-1e-3)).^7/factorial(7)).*((t-tk)>=1e-3);
dt=1e-6;
prc1  = repmat(prc1,1,2);
prc2 = repmat(prc2,1,2);
PIF_tk1=t(PIF_tk1);
PIF_tk2=t(PIF_tk2);

q1 = zeros(1,length(PIF_tk1)-1);
for i = 1:length(PIF_tk1) - 1
    t_integral = PIF_tk1(i):dt:PIF_tk1(i+1);
  
    tk1_pick = PIF_tk1(PIF_tk1 <= PIF_tk1(i));
    temp11 = 0;
    for j = 1:length(tk1_pick)
        temp11 = temp11 + trapz(t_integral, h11(t_integral, tk1_pick(j)).* prc1([1:length(t_integral)]));
    end
    
    tk2_pick = PIF_tk2(PIF_tk2 < PIF_tk1(i+1));
    temp21 = 0;
    for j = 1:length(tk2_pick)
        temp21 = temp21 + trapz(t_integral, h21(t_integral, tk2_pick(j)).* prc1([1:length(t_integral)]));
    end
    q1(i) = thr1 - bias1 * (PIF_tk1(i+1) - PIF_tk1(i)) - temp11 - temp21;
end

q2 = zeros(1,length(PIF_tk2)-1);
for i = 1:length(PIF_tk2) - 1
    t_integral = PIF_tk2(i):dt:PIF_tk2(i+1);
    
    tk1_pick = PIF_tk1(PIF_tk1 < PIF_tk2(i+1));
    temp12 = 0;
    for j = 1:length(tk1_pick)
        temp12 = temp12 + trapz(t_integral, h12(t_integral, tk1_pick(j)).* prc2([1:length(t_integral)]));
    end
    
    tk2_pick = PIF_tk2(PIF_tk2 <= PIF_tk2(i));
    temp22 = 0;
    for j = 1:length(tk2_pick)
        temp22 = temp22 + trapz(t_integral, h22(t_integral, tk2_pick(j)).* prc2([1:length(t_integral)]));
    end

    q2(i) = thr2 - bias2 * (PIF_tk2(i+1) - PIF_tk2(i)) - temp12 - temp22;
end


end