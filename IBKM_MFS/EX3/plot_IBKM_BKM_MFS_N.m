load('NBH90_IBKM_BKM_MFS.mat')
m_bkm=NBH90_BKM;
m_hg=NBH90_IBKM_MFS;
r=m_hg(2:end,1);
nb=m_bkm(1,2:end);
%%BKM
semilogy(nb,m_bkm(2,2:end),'k-*','LineWidth',2);
hold on
for i_g=14   
semilogy(nb,m_hg(i_g,2:end),'r-*','LineWidth',1);
hold on
end
for i_g=96    %%%MFS
semilogy(nb,m_hg(i_g,2:end),'b-x','LineWidth',1.7);
hold on
end

xlabel('NB1');
ylabel('\epsilon_{{R}}')
legend({'BKM','IBKM(R=9.5)','MFS(r0=6.5, r1=0.08, r2=0.02)'},'FontSize',12,'Fontname', 'Times New Roman')
set(gca,'FontSize',12,'Fontname', 'Times New Roman');
box off
set(gca,'LineWidth',1.5)
axis([800 1200 1E-15 1E+5])