
% ���ɹ�һ��ƽ�����е���ʸ��
% �����֤������֤����ȷ
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981
function [a] = planar_streering_vector(Nv,Nh,lambda,d,theta,phi)
% Nv����ֱ��������
% Nh��ˮƽ��������
% lambda���ز�����
% d����Ԫ���
% theta���춥��
% phi����λ��
% a��ƽ�����еĹ�һ������ʸ������������NvNh��1��

% ���ɴ�ֱ����ʸ��
nv= -Nv/2:1:Nv/2-1;
av= 1/sqrt(Nv)*exp(1i*nv.*(2*pi*d/lambda*cos(theta)));
av= av.';

% ���ɴ�ֱ����ʸ��
nh= -Nh/2:1:Nh/2-1;
ah= 1/sqrt(Nh)*exp(1i*nh.*(2*pi*d/lambda*sin(theta)*cos(phi)));
ah= ah.';

% ����ƽ�����е���ʸ��
a= kron(av,ah); 

end
