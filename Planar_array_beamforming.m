
% ƽ������������β���ͼ����������+��λ�ף��ĺ�ģģ�Ⲩ������
% �����֤������֤����ȷ
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981
clc
clear
close all

Nv=8;  % ƽ�����д�ֱ��������
Nh=128; % ƽ������ˮƽ��������
fc=60e9; % �ز�Ƶ��60GHz
c=3e8; % ����
lambda=c/fc; % �ز�����
d=lambda/2; % ������Ԫ���

% ����ͼ������ֵ
theta= pi/180*120; % ������
A0=[4 7 8 5 3 9 6 8];
phi= pi/180*[40 50 70 90 100 120 130 150]; % ��λ��
phi_min= phi- pi/180*1;
phi_max= phi+ pi/180*1;
fai_nv= 2*pi*[0.2 0.4 0.7 0.5 0.9 0.6 0.5 0.3]; % Ƶ����λ

% ����ó���������ʸ��
f= zeros(Nv*Nh,1); % ��������ʸ����ʼ��
mu= 0; % ��һ��ϵ����ʼ��
nv=-Nv/2:1:Nv/2-1;
fai0= fai_nv+ 2*pi*d/lambda*cos(theta)*nv;
for i=1:Nv
    [Fnv,fnv,~]= ULA_beamforming_vector(Nh,A0(i),fai0(i),theta,phi_min(i),phi_max(i),lambda,d); % ���ú���
    f(1+(i-1)*Nh:i*Nh)= Fnv*fnv;
    mu= mu+ Fnv^2;
end
mu= sqrt(mu);
f= 1/mu*f;
    
% ������֤����ͼ��
N_sim= 10240;
g= zeros(1,N_sim); % ����ͼ����ʼ��
phi_sim= 0:pi/N_sim:pi-pi/N_sim;
for i=1:N_sim
    [a]= planar_streering_vector(Nv,Nh,lambda,d,theta,phi_sim(i)); 
    g(i)= a'*f;
end

Loc=[];
% ��ͼ
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
xlabel('\phi(��)');
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
xlabel('\phi(��)');
ylabel('Phase Spectrum (rad)');
legend('simulated','theoretical')
axis([0 180 -pi pi]);
grid on
    