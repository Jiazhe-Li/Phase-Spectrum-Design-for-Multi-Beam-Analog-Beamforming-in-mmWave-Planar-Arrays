
% ���ɼ�Ȩģ�Ⲩ������ʸ�������Ӧ������β���ͼ��
% �����֤������֤����ȷ
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981
function [F,f,phase_spec] = ULA_beamforming_vector(Nh,A0,fai0,theta,phi_min,phi_max,lambda,d)
% Nh��ƽ������ˮƽ��������
% A0�����β���ͼ���ķ�ֵ
% fai0�����β���ͼ���ĺ㶨��λ
% theta�����β���ͼ�����춥��
% phi_min�����β���ͼ������С��λ��
% phi_max�����β���ͼ�������λ��
% lambda���ز�����
% d����Ԫ���
% F��ģ�Ⲩ������ʸ���ļ�Ȩϵ��
% f����һ��ģ�Ⲩ������ʸ������������Nh��1��

% ���ֽ�Ƶ�����䷶Χ
omega_min= 2*pi*d/lambda*sin(theta)*cos(phi_max);
omega_max= 2*pi*d/lambda*sin(theta)*cos(phi_min);

% ģ�Ⲩ������ʸ���ļ�Ȩϵ��
F= sqrt(A0^2/(2*pi)*(omega_max-omega_min));
C= -A0^2*Nh/(4*pi*F^2)*(omega_min+omega_max);

% ULA���е���λ
nh= -Nh/2:1:Nh/2-1;
fai= pi*F^2/(A0^2*Nh)*(nh.^2-2*C*nh)+ fai0;

% ULA���еĹ�һ����������ʸ��
f= 1/sqrt(Nh)*exp(1i*fai);
f= f.';

% ���㲨��ͼ������λ��
fai_d= pi*F^2/(A0^2*Nh)*(2*nh-2*C); % fai�ĵ���
phase_spec= angle( exp( 1i*(fai- fai_d.*nh)+ 1i*0.7  ) );

end
