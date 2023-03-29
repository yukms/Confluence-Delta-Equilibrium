%% input parameter
clear
H= 6; % 낙동강 깊이
L=360; % 낙동강 너비
S=230; % 현재 Delta progradation (평형상태로 가정)
B=320;
qs = 43563; % t/ year 
qs = qs * (1000^2) / 2.65 / (B * 100) / 10000 ; % m2/yr
qs = qs / 31536000; % m2 / sec
Qw = 200; % m3/s
% alpha = (L-S)* qs/ Qw;    % [  ] 무차원
alpha =4.54E-06;
St = 0.001;  % topset slope
Sf = 0.01 ; % foreset slope
dt = 60*60*24;

%% alpha model 삼각주의 시간에 따른 변화와 equlibrium 상태로의 수렴
% Main equation
% St * s(i)^2 + (H- s(i-1)*St) * s(i)   - s(i-1)*H - qs*dt + Edot*dt*H
% Edot = alpha * Qw / H / (L-s(i))
tend=350000; % day 
duration=30; % duration of event
extreme=20;
Qw=ones(tend,1)*Qw;
qs=ones(tend,1)*qs;
s=zeros(tend,1);
Edot=zeros(tend,1);
Area=zeros(tend,1);
s_i = 0;   
Qw(50000:50000+duration)=Qw(1)*extreme;   % Qw 증가
Qw(100000:100000+duration*5)=Qw(1)/extreme;  % Qw 감소   % 감소 duration은 좀 길게~
qs(150000:150000+duration)=qs(1)*extreme;  % qs 증가
qs(200000:200000+duration*5)=qs(1)/extreme;  % qw 감소

Qw(250000:250000+duration)=Qw(1)*extreme; % 모두증가
qs(250000:250000+duration)=qs(1)*extreme;
Qw(300000:300000+duration*5)=Qw(1)/extreme; % 모두감소
qs(300000:300000+duration*5)=qs(1)/extreme;

for i = 1 : tend
    if i == 1
    Area(i) = H* (360-s_i); 
    Edot(i) = alpha * (Qw(i)^2) / (Area(i)^2);
    s(i) = ((s_i*St-H) + sqrt (  (H-s_i*St)^2 - 4 *St *(-s_i*H-qs(i)*dt+Edot(i)*dt*H)))/(2*St);
    else
    Area(i) = H* (360-s(i-1));
    Edot(i) = alpha * (Qw(i)^2) / (Area(i)^2);
    s(i) = ((s(i-1)*St-H) + sqrt (  (H-s(i-1)*St)^2 - 4 *St *(-s(i-1)*H-qs(i)*dt+Edot(i)*dt*H)))/(2*St);
    end
    
end

figure(9)   %% figure9
subplot(5,1,1)
semilogy(Qw,'k')
xlim([0 tend])
ylabel('$Q_w [m^3/s]$','Interpreter','latex','FontSize',12)

subplot(5,1,2)
semilogy(qs,'k')
ylabel('$q_s [m^2/s]$','Interpreter','latex','FontSize',12)
xlim([0 tend])


subplot(5,1,3)
plot(s,'k')
ylabel('$S [m]$','Interpreter','latex','FontSize',12)

subplot(5,1,4)
plot(Area,'k')
ylabel('Area $[m^2]$','Interpreter','latex','FontSize',12)
subplot(5,1,5)
semilogy(Edot,'k')
ylabel('$\dot{E}$','Interpreter','latex','FontSize',12)
xlabel('day')

%%
filename = 'gam_qs_gumi_Qw.csv';
T1 = readtable(filename);
xtime=T1.x__(1:2557);

T=readmatrix('gam_qs_gumi_Qw.csv');
qs_gam=T(:,7);
TF= isnan(qs_gam);
qs_gam(TF)=0;
nantomean=mean(qs_gam);
qs_gam(TF)=nantomean;
Qw_bo=T(:,8);
TF2=isnan(Qw_bo);
Qw_bo(TF2)=0;
nantomean2=mean(Qw_bo);
Qw_bo(TF2)=nantomean2;
H= 6; % 낙동강 깊이
L=360; % 낙동강 너비
S=230; % 현재 Delta progradation (평형상태로 가정)
B=320;
tend=length(qs_gam); % 2557 day

startDate = datenum('01-01-2015');
endDate = datenum('12-31-2021');


s=zeros(tend,1);
Edot=zeros(tend,1);
Area=zeros(tend,1);
s_i = 210;   

for i = 1 : tend
   
    if i == 1
    Area(i) = H* (360-s_i); 
    Edot(i) = alpha * (Qw_bo(i)^2) / (Area(i)^2);
    s(i) = ((s_i*St-H) + sqrt (  (H-s_i*St)^2 - 4 *St *(-s_i*H-qs_gam(i)*dt+Edot(i)*dt*H)))/(2*St);
    else
    Area(i) = H* (360-s(i-1));
    Edot(i) = alpha * (Qw_bo(i)^2) / (Area(i)^2);
    s(i) = ((s(i-1)*St-H) + sqrt (  (H-s(i-1)*St)^2 - 4 *St *(-s(i-1)*H-qs_gam(i)*dt+Edot(i)*dt*H)))/(2*St);
    end
    
end

figure(10)  %% figure 10
subplot(5,1,1)
semilogy(xtime,Qw_bo,'k')
ylabel('$Q_w [m^3/s]$','Interpreter','latex','FontSize',12)

subplot(5,1,2)
semilogy(xtime,qs_gam,'k')
ylabel('$q_s [m^2/s]$','Interpreter','latex','FontSize',12)


subplot(5,1,3)
plot(xtime,s,'k')
ylabel('$S [m]$','Interpreter','latex','FontSize',12)
xlabel('day')
hold on
plot( [17 94 216 729 1049 1189 1250 1254 1381 1810 ],[211 215 233 260  204 205 222 222 186 210],'r+','MarkerSize',10)
yline(226.2401,'r--','linewidth',1)
legend('Model Prediction','Mearsured S','Model Predicted equilibrium S')

subplot(5,1,4)
plot(xtime,Area,'k')
ylabel('Area $[m^2]$','Interpreter','latex','FontSize',12)
subplot(5,1,5)
semilogy(xtime,Edot,'k')
ylabel('$\dot{E}$','Interpreter','latex','FontSize',12)


