
% 生成归一化平面阵列导向矢量
% 结果验证：已验证，正确
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981
function [a] = planar_streering_vector(Nv,Nh,lambda,d,theta,phi)
% Nv：垂直天线数量
% Nh：水平天线数量
% lambda：载波波长
% d：阵元间距
% theta：天顶角
% phi：方位角
% a：平面阵列的归一化导向矢量，列向量，NvNh行1列

% 生成垂直导向矢量
nv= -Nv/2:1:Nv/2-1;
av= 1/sqrt(Nv)*exp(1i*nv.*(2*pi*d/lambda*cos(theta)));
av= av.';

% 生成垂直导向矢量
nh= -Nh/2:1:Nh/2-1;
ah= 1/sqrt(Nh)*exp(1i*nh.*(2*pi*d/lambda*sin(theta)*cos(phi)));
ah= ah.';

% 生成平面阵列导向矢量
a= kron(av,ah); 

end
