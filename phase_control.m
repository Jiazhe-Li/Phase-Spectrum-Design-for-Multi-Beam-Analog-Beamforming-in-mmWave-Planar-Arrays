clc
clear
close all

% ƽ������������β���ͼ����������+��λ�ף��ĺ�ģģ�Ⲩ������
% �����֤������֤����ȷ
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981
clc
clear
close all

Nv=1;  % ƽ�����д�ֱ��������
Nh=4096; % ƽ������ˮƽ��������
fc=60e9; % �ز�Ƶ��60GHz
c=3e8; % ����
lambda=c/fc; % �ز�����
d=lambda/2; % ������Ԫ���

% ����ͼ������ֵ
theta= pi/180*90; % ������
A0=1;
phi= pi/180*120; % ��λ��
phi_min= phi- pi/180*20;
phi_max= phi+ pi/180*20;
fai_nv= 2*pi*0.2; % Ƶ����λ

% ����ó���������ʸ��
f= zeros(Nv*Nh,1); % ��������ʸ����ʼ��
mu= 0; % ��һ��ϵ����ʼ��
nv=-Nv/2:1:Nv/2-1;
fai0= fai_nv+ 2*pi*d/lambda*cos(theta)*nv;
for i=1:Nv
    [Fnv,fnv,spec]= ULA_beamforming_vector(Nh,A0(i),fai0(i),theta,phi_min(i),phi_max(i),lambda,d); % ���ú���
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
plot(180/pi*phi_sim,angle(g_));
hold on
for i=1:length(fai_nv)
%     plot(180/pi*phi(i),spec,'ro');
    wmin= 2*pi*d/lambda*sin(theta)*cos(phi_max);
    wmax= 2*pi*d/lambda*sin(theta)*cos(phi_min);
    plot(acos( lambda/(2*pi*d)*(wmin:(wmax-wmin)/Nh:wmax-(wmax-wmin)/Nh) ).*(180/pi),spec,'r'); % ���ھ��β���ͼ����w��n�����Թ�ϵ�����Կ���ֱ�Ӿ��ȼ����w�������꣬������Ҫ�󷴺���
end
xlabel('\phi(��)');
ylabel('Phase Spectrum (rad)');
legend('simulated','theoretical')
axis([110 130 -pi pi]);
grid on
    