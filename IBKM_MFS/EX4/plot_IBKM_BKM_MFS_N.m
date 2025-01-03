load('NBH60_IBKM_BKM_MFS.mat')
m_bkm=NBH60_BKM;
m_hg=NBH60_IBKM_MFS;
r=m_hg(2:end,1);
nb=m_bkm(1,2:end);
%%BKM
semilogy(nb,m_bkm(2,2:end),'k-*','LineWidth',2);
hold on
for i_g=52     %%%R=12
semilogy(nb,m_hg(i_g,2:end),'r-*','LineWidth',1);
 hold on
end
for i_g=50   %%%MFS
semilogy(nb,m_hg(i_g,2:end),'b-x','LineWidth',1.7);
 hold on
end

xlabel('NB1');
ylabel('\epsilon_{{R}}')
legend({'BKM','IBKM-Cloud(R=12)','MFS(r0=3.5, r1=0.18, r2=0.12, r3=0.05)'},'FontSize',12,'Fontname', 'Times New Roman')
set(gca,'FontSize',12,'Fontname', 'Times New Roman');
box off
set(gca,'LineWidth',1.5)
axis([800 1200 1E-15 1E+5])