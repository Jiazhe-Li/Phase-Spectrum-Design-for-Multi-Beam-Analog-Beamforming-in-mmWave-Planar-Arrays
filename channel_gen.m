
% 产生毫米波稀疏信道
clc
clear
close all
% 结果验证：已验证，正确（包括各个簇、各条子径的功率验证）
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981

Nv=8;  % 平面阵列垂直天线数量
Nh=256; % 平面阵列水平天线数量
fc=60e9; % 载波频率60GHz
c=3e8; % 光速
lambda=c/fc; % 载波波长
d=lambda/2; % 天线阵元间距

% 生成路径簇的数量L
L=min(max(1,random('Poisson',1.6)),4); % 服从泊松分布

% 生成路径簇的功率
r_tao=2.8;
sita=4.0;
Z=randn(1,L)*sita;
U=rand(1,L);
gamma=(U.^(r_tao-1)).*10.^(-0.1*Z);
gamma=gamma./sum(gamma); % 功率归一化

% 生成路径簇的中心角度
phi_min= pi/6+ (2/3*pi)/L*((1:L)-1);
phi_max= pi/6+ (2/3*pi)/L*(1:L);
phi= rand(1,L).*(phi_max-phi_min)+ phi_min; % 每个路径簇的角度，在(phi_min,phi_max)内服从均匀分布

% 生成路径簇中每条子径的角度
sigma_rms= pi/180*10; % 假设每条路径簇的RMS角度扩展为10°
P=20; % 假设每条路径簇中含有P条子径
Phi= zeros(L,P); % 所有子径角度初始化
for i=1:L
    Phi(i,:)= mod( phi(i)+sigma_rms*randn(1,P), pi ); % 子径角度服从缠绕正态分布，定义域为(0,pi)
end

% 生成每条子径的导向矢量
theta= pi/180*120; % 天顶角固定为120°
A= cell(1,L);
for i=1:L
    for p=1:P
        A{1,i}(:,p)=planar_streering_vector(Nv,Nh,lambda,d,theta,Phi(i,p)); % 信道的导向矢量（已归一化，列矢量）
    end
end

% 生成每条子径的复增益
Alpha= zeros(L,P);
for i=1:L
    for p=1:P
        Alpha(i,p)= (randn(1)+1i*randn(1))/sqrt(2)*sqrt(gamma(i)/P);
    end
end

% 合并成多径信道
h=zeros(Nv*Nh,1);
for i=1:L
    for p=1:P
        h= h+ Alpha(i,p)*A{1,i}(:,p);
    end
end
h= sqrt(Nv*Nh)*h; %  功率归一化，即满足 E(h*h')=I

% 仿真信道的波束图案
N_sim= 10240;
g= zeros(1,N_sim); % 波束图案初始化
phi_sim= 0:pi/N_sim:pi-pi/N_sim;
for i=1:N_sim
    [a]= planar_streering_vector(Nv,Nh,lambda,d,theta,phi_sim(i)); 
    g(i)= h'*a;
end
  
% 画图
figure(1)
subplot(2,1,1)
plot(180/pi*phi_sim,abs(g));
xlabel('方位角(°)');
ylabel('信道增益');
grid on
axis tight
subplot(2,1,2)
g_=g;
g_(abs(g_)<3)=0;
plot(180/pi*phi_sim,angle(g_)./(2*pi));
xlabel('方位角(°)');
ylabel('相位');
grid on
axis tight
