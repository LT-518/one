load('ibkm_mfs.mat')
m_bkm=ibkm;
m_hg=mfs;
nb=m_bkm(1,2:end);
for i_g=68     %ibkm
semilogy(nb,m_bkm(i_g,2:end),'r-*','LineWidth',1);
hold on
end
%%mfs
semilogy(nb,m_hg(2,2:end),'b-*','LineWidth',2);
hold on

xlabel('NB1');
ylabel('\epsilon_{{R}}')
legend({'IBKM(R=9.6)','MFS(r0=4.5, r1=0.008)'},'FontSize',12,'Fontname', 'Times New Roman')
set(gca,'FontSize',12,'Fontname', 'Times New Roman');
box off
set(gca,'LineWidth',1.5)
axis([100 400 3.5E-6 1E-3])