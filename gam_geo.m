%% model run
for i = 1 : 1 : tend
   
    
    if i ==1 
        S(i) = (-H+sqrt(H^2-2*St*(-H*S_initial-0.5*St*S_initial^2 -qs_gam)))/St;
        progradationrate(i)=S(i);

    else
        S(i) = (-H+sqrt(H^2-2*St*(-H*S(i-1)-0.5*St*S(i-1)^2 -qs_gam)))/St;
        R(i,1) = B-S(i,1);
        A(i,1) = R(i)*H;
        v(i,1) = Qw_nakdong /A(i);
                       progradationrate(i,1)=S(i)-S(i-1);

            if v(i,1) > v_critical
            S(i,1)=S(i,1)-v(i,1)*erosioncoef;
            erosion(i,1)=erosioncoef*v(i,1);
            end
                       progradationrate2(i,1)=S(i)-S(i-1);

            
    end

    


end
f=figure;
f.Position=[1200 0 600 1000]
subplot(3,1,1)
hold on
plot(S)
yline(350)
ylabel('shoreline [m]')
subplot(3,1,2)
plot(v)
ylabel('velocity [m/s]')
hold on
yline(v_critical)
subplot(3,1,3)
plot(erosion,'r')
hold on
plot(progradationrate,'b--')
% plot(progradationrate2,'g:');

% plot(progradationrate-erosion)
xlim([2 tend])
legend('erosion rate','progradation rate','progradtaion rate +erosion','Location','southwest')