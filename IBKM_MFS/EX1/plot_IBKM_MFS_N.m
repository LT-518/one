load('ibkm_mfs.mat')
m_bkm=ibkm_mfs;
nb=m_bkm(1,2:end);

for i_g=3    %mfs
semilogy(nb,m_bkm(i_g,2:end),'b-*','LineWidth',1);
hold on
end
for i_g=2     %ibkm
semilogy(nb,m_bkm(i_g,2:end),'r-*','LineWidth',1);
hold on
end
xlabel('NB1');
ylabel('\epsilon_{{R}}')
legend({'IBKM(R=0.06)','MFS(r0=0.09, r1=0.06)'},'FontSize',12,'Fontname', 'Times New Roman')
set(gca,'FontSize',12,'Fontname', 'Times New Roman');
box off
set(gca,'LineWidth',1.5)
axis([40 200 1E-3 4E-1])