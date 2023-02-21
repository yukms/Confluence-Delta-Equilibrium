%% input parameter

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

%% alpha sensitivity analysis : S일정할 떄 qs Qw에 따른 alpha 변화
qs_range = [1000 : 100 : 90000-9000];
qs_range = qs_range * (1000^2) / 2.65 / (320 * 100) / 10000 ; % m2/yr
qs_range = qs_range / 31536000; % m2 / sec
qs_range=qs_range';
qs_range_inv= 1./qs_range;  % 행렬의 나누기? 를 그냥 곱셈으로 사용하기 위해 inverse
Qw_range = [50:0.5:450];
Qw_range_inv = 1./Qw_range;
% alpha_range = (L-S)*qs_range*Qw_range_inv; % U1
alpha_range = (L-S)^2*qs_range*6*(Qw_range_inv).^2; % U2

alpha_range(alpha_range<1E-07)=nan;

% contourf figure 
figure(1)
contourf(Qw_range,qs_range,alpha_range);
colorbar;
c = colorbar;
c.Label.String = '\alpha  erosion coefficient';
xlabel('Q_w [m^2/sec]')
ylabel('q_s [m^2/sec]')
hold on
plot(Qw,qs,'rp' ,'MarkerFaceColor','red','MarkerSize',15)

% pcolor interp figure
figure(2) 
pcolor(Qw_range,qs_range,alpha_range)
hold on
shading interp
contourf(alpha_range) % not working?
colorbar;
c = colorbar;
c.Label.String = '\alpha  erosion coefficient';
xlabel('Q_w [m^2/sec]')
ylabel('q_s [m^2/sec]')
hold on
plot(Qw,qs,'rp' ,'MarkerFaceColor','red','MarkerSize',15)
%% S sensitivity analysis : alpha 일정할 때 qs Qw에 따른 equil S의 변화
% S_range = L- alpha * qs_range_inv*Qw_range ;  % U1
S_range = L -  sqrt(alpha/6./qs_range)*Qw_range; % U2
S_range(S_range<0)=nan;
% contourf figure
figure(3)
contourf(Qw_range,qs_range,S_range);
colorbar;
c2=colorbar;
c2.Label.String = 'equilibrium shoreline location [m]';
xlabel('Q_w [m^2/sec]')
ylabel('q_s [m^2/sec]')
hold on
% plot(Qw,qs,'rp' ,'MarkerFaceColor','red','MarkerSize',15)

% smoothing figure
figure(4)
pcolor(Qw_range,qs_range,S_range)
hold on
shading interp
contourf(S_range) % not working?
c2=colorbar;
c2.Label.String = 'equilibrium shoreline location [m]';
xlabel('Q_w [m^2/sec]')
ylabel('q_s [m^2/sec]')
hold on
plot(Qw,qs,'rp' ,'MarkerFaceColor','red','MarkerSize',15)

% figures for alpha, S sensitivity (smooth) 
figure(5)
subplot(2,1,1)
pcolor(Qw_range,qs_range,alpha_range)
hold on
shading interp
contourf(alpha_range) % not working?
colorbar;
c = colorbar;
c.Label.String = '\alpha  erosion coefficient';
xlabel('Q_w [m^2/sec]')
ylabel('q_s [m^2/sec]')
hold on
plot(Qw,qs,'rp' ,'MarkerFaceColor','red','MarkerSize',15)
subplot(2,1,2)
pcolor(Qw_range,qs_range,S_range)
hold on
shading interp
contourf(S_range) % not working?
c2=colorbar;
c2.Label.String = 'equilibrium shoreline location [m]';
xlabel('Q_w [m^2/sec]')
ylabel('q_s [m^2/sec]')
hold on
plot(Qw,qs,'rp' ,'MarkerFaceColor','red','MarkerSize',15)


    


%% alpha 모든 합류부 일정하다고 가정하였을 때 다른 델타들의 거동을 예상
figure
set( gcf, 'Position', [0 0 800 600] ) ;
S_range(S_range<0)=nan;
S_range(S_range>L)=nan;
contourf(Qw_range,qs_range,S_range/L);
colorbar;
c2=colorbar;
c2.Label.String = 'equilibrium shoreline location [m]';
xlabel('Q_w [m^2/sec]','FontSize',13)
ylabel('q_s [m^2/sec]','FontSize',13)
hold on
plot(Qw,qs,'rp' ,'MarkerFaceColor','red','MarkerSize',10)

plot(73.09,2.18E-06,'bp','MarkerSize',10,'MarkerFaceColor','b')
plot(150.7828,1.53E-06,'o','MarkerSize',10,'MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840])
plot(278.0292,2.73E-07,'c<','MarkerSize',10,'MarkerFaceColor','c')
plot(289.2921,2.18E-06,'y^','MarkerSize',10,'MarkerFaceColor','y')
plot(337.8056,5.49E-07,'md','MarkerSize',10,'MarkerFaceColor','m')
legend('sensitivity range','감천','내성천','병성천','회천','황강','남강','FontSize',10)
    

%% 추세선
uu(:,1)=[0.092769129,0.011938381, 0.063310173, 0.023858098, 0.055317548, 0.025875878];
uu(:,2)=[3.64E-07, 2.54E-07, 2.71E-07, 4.55E-08,3.63E-07,9.14E-08];
i = 0 : 0.01 : 0.1;
iy= alpha * i;
figure
plot(uu(:,1),uu(:,2),'ok','MarkerSize',5,'MarkerFaceColor','k')
hold on
plot(i,iy,'k')
legend('Nakdong RIver Confluences','Linear Regression (y=4.54E-06x)','Location','northwest','FontSize',12)
xlabel('U^2','FontSize',15)
ylabel('$\dot{E}$','Interpreter','latex','FontSize',15)

%%
bardata = [0.692471892 0.714285714; 0.753057562   0.465116279 ; 0.628310955  0.638888889 ; 0  0 ; 0.544829257  0.453333333; 0 0]
bar(bardata)
legend('Model Prediction','Observed')
ylabel('S/L')






    