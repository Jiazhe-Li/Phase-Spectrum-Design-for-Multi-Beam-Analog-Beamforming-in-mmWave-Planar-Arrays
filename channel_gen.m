
% �������ײ�ϡ���ŵ�
clc
clear
close all
% �����֤������֤����ȷ�����������ء������Ӿ��Ĺ�����֤��
% Author: Jiazhe Li, 2025, ORCID 0000-0002-1042-9981

Nv=8;  % ƽ�����д�ֱ��������
Nh=256; % ƽ������ˮƽ��������
fc=60e9; % �ز�Ƶ��60GHz
c=3e8; % ����
lambda=c/fc; % �ز�����
d=lambda/2; % ������Ԫ���

% ����·���ص�����L
L=min(max(1,random('Poisson',1.6)),4); % ���Ӳ��ɷֲ�

% ����·���صĹ���
r_tao=2.8;
sita=4.0;
Z=randn(1,L)*sita;
U=rand(1,L);
gamma=(U.^(r_tao-1)).*10.^(-0.1*Z);
gamma=gamma./sum(gamma); % ���ʹ�һ��

% ����·���ص����ĽǶ�
phi_min= pi/6+ (2/3*pi)/L*((1:L)-1);
phi_max= pi/6+ (2/3*pi)/L*(1:L);
phi= rand(1,L).*(phi_max-phi_min)+ phi_min; % ÿ��·���صĽǶȣ���(phi_min,phi_max)�ڷ��Ӿ��ȷֲ�

% ����·������ÿ���Ӿ��ĽǶ�
sigma_rms= pi/180*10; % ����ÿ��·���ص�RMS�Ƕ���չΪ10��
P=20; % ����ÿ��·�����к���P���Ӿ�
Phi= zeros(L,P); % �����Ӿ��Ƕȳ�ʼ��
for i=1:L
    Phi(i,:)= mod( phi(i)+sigma_rms*randn(1,P), pi ); % �Ӿ��Ƕȷ��Ӳ�����̬�ֲ���������Ϊ(0,pi)
end

% ����ÿ���Ӿ��ĵ���ʸ��
theta= pi/180*120; % �춥�ǹ̶�Ϊ120��
A= cell(1,L);
for i=1:L
    for p=1:P
        A{1,i}(:,p)=planar_streering_vector(Nv,Nh,lambda,d,theta,Phi(i,p)); % �ŵ��ĵ���ʸ�����ѹ�һ������ʸ����
    end
end

% ����ÿ���Ӿ��ĸ�����
Alpha= zeros(L,P);
for i=1:L
    for p=1:P
        Alpha(i,p)= (randn(1)+1i*randn(1))/sqrt(2)*sqrt(gamma(i)/P);
    end
end

% �ϲ��ɶྶ�ŵ�
h=zeros(Nv*Nh,1);
for i=1:L
    for p=1:P
        h= h+ Alpha(i,p)*A{1,i}(:,p);
    end
end
h= sqrt(Nv*Nh)*h; %  ���ʹ�һ���������� E(h*h')=I

% �����ŵ��Ĳ���ͼ��
N_sim= 10240;
g= zeros(1,N_sim); % ����ͼ����ʼ��
phi_sim= 0:pi/N_sim:pi-pi/N_sim;
for i=1:N_sim
    [a]= planar_streering_vector(Nv,Nh,lambda,d,theta,phi_sim(i)); 
    g(i)= h'*a;
end
  
% ��ͼ
figure(1)
subplot(2,1,1)
plot(180/pi*phi_sim,abs(g));
xlabel('��λ��(��)');
ylabel('�ŵ�����');
grid on
axis tight
subplot(2,1,2)
g_=g;
g_(abs(g_)<3)=0;
plot(180/pi*phi_sim,angle(g_)./(2*pi));
xlabel('��λ��(��)');
ylabel('��λ');
grid on
axis tight
