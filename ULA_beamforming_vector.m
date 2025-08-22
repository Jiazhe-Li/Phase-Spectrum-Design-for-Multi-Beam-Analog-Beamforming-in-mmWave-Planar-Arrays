
% 生成加权模拟波束赋形矢量，其对应任意矩形波束图案
% 结果验证：已验证，正确
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981
function [F,f,phase_spec] = ULA_beamforming_vector(Nh,A0,fai0,theta,phi_min,phi_max,lambda,d)
% Nh：平面阵列水平天线数量
% A0：矩形波束图案的幅值
% fai0：矩形波束图案的恒定相位
% theta：矩形波束图案的天顶角
% phi_min：矩形波束图案的最小方位角
% phi_max：矩形波束图案的最大方位角
% lambda：载波波长
% d：阵元间距
% F：模拟波束赋形矢量的加权系数
% f：归一化模拟波束赋形矢量，列向量，Nh行1列

% 数字角频率区间范围
omega_min= 2*pi*d/lambda*sin(theta)*cos(phi_max);
omega_max= 2*pi*d/lambda*sin(theta)*cos(phi_min);

% 模拟波束赋形矢量的加权系数
F= sqrt(A0^2/(2*pi)*(omega_max-omega_min));
C= -A0^2*Nh/(4*pi*F^2)*(omega_min+omega_max);

% ULA阵列的相位
nh= -Nh/2:1:Nh/2-1;
fai= pi*F^2/(A0^2*Nh)*(nh.^2-2*C*nh)+ fai0;

% ULA阵列的归一化波束赋形矢量
f= 1/sqrt(Nh)*exp(1i*fai);
f= f.';

% 推算波束图案的相位谱
fai_d= pi*F^2/(A0^2*Nh)*(2*nh-2*C); % fai的导数
phase_spec= angle( exp( 1i*(fai- fai_d.*nh)+ 1i*0.7  ) );

end
