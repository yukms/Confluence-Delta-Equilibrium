clear
close all
%% confluence delta model


%% Input parameters
Qw_nakdong = 10;
qs_gam = 10; % gam river로 부터 공급되는 qs ; Qs/B (1~11)
erosioncoef=3000;
v_critical=0.01; % v > v_critical  >>  washover (delta retreat)
S_initial = 10; % initial delta progradation; in meters


tend = 300;
H = 6; % water depth in meters
B = 350; % width of nakdong river; in meters (normal stage ; interflood)
v = zeros( tend,1);
R = zeros( tend,1);
S = zeros( tend,1);
A = zeros ( tend,1);
St=0.001;
Sf=0.01;
%% 추후에 어떤 shear stress,, wash over equation 등 적용 가능

gam_geo