NB1 = [800, 900, 1000, 1100,1200, 800, 900, 1000, 1100,1200];
NBH = [90, 90, 90, 90, 90,120, 120, 120, 120,120];
CPU_ECN = [4.4813, 3.6504, 5.0285, 11.2347, 15.6342,7.2234 ,8.2619 ,11.6717, 15.1133,19.2940];
CPU_EECN = [2.7576,3.3764 ,4.0501, 5.4349, 	6.8236,3.2727, 4.2386 ,	5.3291,7.0107,8.6946 ];

NB1_90 = NB1(NBH == 90);
CPU_ECN_90 = CPU_ECN(NBH == 90);
CPU_EECN_90 = CPU_EECN(NBH == 90);

NB1_120 = NB1(NBH == 120);
CPU_ECN_120 = CPU_ECN(NBH == 120);
CPU_EECN_120 = CPU_EECN(NBH == 120);

figure('Color', 'w');
%plot NBH=90
yyaxis left;
plot(NB1_90, CPU_ECN_90, '-o', 'DisplayName', 'ECN (NBH=90)');
hold on;
plot(NB1_90, CPU_EECN_90, '-x', 'DisplayName', 'EECN (NBH=90)');
ylabel('CPU Run Time (seconds, NBH=90)');

% plot NBH=120
yyaxis right;
plot(NB1_120, CPU_ECN_120, '--o', 'DisplayName', 'ECN (NBH=120)');
plot(NB1_120, CPU_EECN_120, '--x', 'DisplayName', 'EECN (NBH=120)');
ylabel('CPU Run Time (seconds, NBH=120)');

legend('show');

title('Comparison of CPU Run Time for ECN and EECN');
xlabel('NB1');
grid off;