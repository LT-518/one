NB1 = [200, 400, 600, 800, 200, 400, 600, 800];
NBH = [40, 40, 40, 40, 80, 80, 80, 80];
CPU_ECN = [0.1948, 0.4603, 0.9350, 1.9115,0.2754,0.5991, 1.2440 ,2.3702 ];
CPU_EECN = [0.1321, 0.2895, 0.5793, 1.0062, 0.1825, 0.3829, 0.7079,1.2724];

NB1_40 = NB1(NBH == 40);
CPU_ECN_40 = CPU_ECN(NBH == 40);
CPU_EECN_40 = CPU_EECN(NBH == 40);

NB1_80 = NB1(NBH == 80);
CPU_ECN_80 = CPU_ECN(NBH == 80);
CPU_EECN_80 = CPU_EECN(NBH == 80);

figure('Color', 'w');
%plot NBH=40
yyaxis left;
plot(NB1_40, CPU_ECN_40, '-o', 'DisplayName', 'ECN (NBH=40)');
hold on;
plot(NB1_40, CPU_EECN_40, '-x', 'DisplayName', 'EECN (NBH=40)');
ylabel('CPU Run Time (seconds, NBH=40)');

% plot NBH=80
yyaxis right;
plot(NB1_80, CPU_ECN_80, '--o', 'DisplayName', 'ECN (NBH=80)');
plot(NB1_80, CPU_EECN_80, '--x', 'DisplayName', 'EECN (NBH=80)');
ylabel('CPU Run Time (seconds, NBH=80)');

legend('show');

title('Comparison of CPU Run Time for ECN and EECN');
xlabel('NB1');
grid off;