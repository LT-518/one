NB1 = [400, 500, 600, 700, 400, 500, 600, 700];
NBH = [60, 60, 60, 60, 80, 80, 80, 80];
CPU_ECN = [1.2543,2.1868,2.9766,4.0108,1.5031,2.3514,3.5012,4.5388];
CPU_EECN = [0.8784,1.2751,1.7401,2.4951,1.0013,1.5571,1.9659,2.6837];

NB1_60 = NB1(NBH == 60);
CPU_ECN_60 = CPU_ECN(NBH == 60);
CPU_EECN_60 = CPU_EECN(NBH == 60);

NB1_80 = NB1(NBH == 80);
CPU_ECN_80 = CPU_ECN(NBH == 80);
CPU_EECN_80 = CPU_EECN(NBH == 80);

figure('Color', 'w');

% plot 60
yyaxis left;
plot(NB1_60, CPU_ECN_60, '-o', 'DisplayName', 'ECN (NBH=60)');
hold on;
plot(NB1_60, CPU_EECN_60, '-x', 'DisplayName', 'EECN (NBH=60)');
ylabel('CPU Run Time (seconds, NBH=60)');

% plot 80
yyaxis right;
plot(NB1_80, CPU_ECN_80, '--o', 'DisplayName', 'ECN (NBH=80)');
plot(NB1_80, CPU_EECN_80, '--x', 'DisplayName', 'EECN (NBH=80)');
ylabel('CPU Run Time (seconds, NBH=80)');

legend('show');

title('Comparison of CPU Run Time for ECN and EECN');
xlabel('NB1');

grid off;