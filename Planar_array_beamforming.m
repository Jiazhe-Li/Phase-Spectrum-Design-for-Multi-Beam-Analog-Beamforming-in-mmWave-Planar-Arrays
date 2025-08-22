
% 平面阵列任意矩形波束图案（幅度谱+相位谱）的恒模模拟波束赋形
% 结果验证：已验证，正确
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981
clc
clear
close all

Nv=8;  % 平面阵列垂直天线数量
Nh=128; % 平面阵列水平天线数量
fc=60e9; % 载波频率60GHz
c=3e8; % 光速
lambda=c/fc; % 载波波长
d=lambda/2; % 天线阵元间距

% 波束图案理论值
theta= pi/180*120; % 俯仰角
A0=[4 7 8 5 3 9 6 8];
phi= pi/180*[40 50 70 90 100 120 130 150]; % 方位角
phi_min= phi- pi/180*1;
phi_max= phi+ pi/180*1;
fai_nv= 2*pi*[0.2 0.4 0.7 0.5 0.9 0.6 0.5 0.3]; % 频谱相位

% 计算得出波束赋形矢量
f= zeros(Nv*Nh,1); % 波束赋形矢量初始化
mu= 0; % 归一化系数初始化
nv=-Nv/2:1:Nv/2-1;
fai0= fai_nv+ 2*pi*d/lambda*cos(theta)*nv;
for i=1:Nv
    [Fnv,fnv,~]= ULA_beamforming_vector(Nh,A0(i),fai0(i),theta,phi_min(i),phi_max(i),lambda,d); % 调用函数
    f(1+(i-1)*Nh:i*Nh)= Fnv*fnv;
    mu= mu+ Fnv^2;
end
mu= sqrt(mu);
f= 1/mu*f;
    
% 仿真验证波束图案
N_sim= 10240;
g= zeros(1,N_sim); % 波束图案初始化
phi_sim= 0:pi/N_sim:pi-pi/N_sim;
for i=1:N_sim
    [a]= planar_streering_vector(Nv,Nh,lambda,d,theta,phi_sim(i)); 
    g(i)= a'*f;
end

Loc=[];
% 画图
figure(1)
subplot(2,1,1)
plot(180/pi*phi_sim,abs(g)*mu*sqrt(Nv*Nh));
hold on
for i=1:length(A0)
    amp= zeros(1,length(g));
    loc= find( phi_sim >= phi_min(i) & phi_sim <= phi_max(i));
    amp(loc)= A0(i);
    plot(180/pi*phi_sim,amp,'r');
    Loc=[Loc loc];
end
xlabel('\phi(°)');
ylabel('Amplitude Spectrum');
legend('simulated','theoretical')
grid on
axis([0 180 -inf inf]);
subplot(2,1,2)
g_=g*mu*sqrt(Nv*Nh);
% g_(abs(g_)<2)=0;
g_2= zeros(1,length(g_));
g_2(Loc)= g_(Loc);
plot(180/pi*phi_sim,angle(g_2));
hold on
for i=1:length(fai_nv)
    plot(180/pi*phi(i),angle( exp(1i*fai_nv(i)) ),'ro');
end
xlabel('\phi(°)');
ylabel('Phase Spectrum (rad)');
legend('simulated','theoretical')
axis([0 180 -pi pi]);
grid on
    