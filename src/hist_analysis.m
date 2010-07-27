close all
a=I_LR-I;
x=-4.5:0.1:4.5;
hist(a(:),x)
axis([min(x) max(x) 0 3000])
set(findobj(gca,'Type','Patch'),'FaceColor','y','EdgeColor','black')
hold on
y=gaussmf(x,[1,0]);
y=y*numel(a)/sum(y);
plot(x,y,'LineWidth',3,'Color','b')
xlabel('\epsilon')
ylabel('f_{LR}(\epsilon)')
legend('Histogram of \epsilon','N(0,1)')
% legend('Histogram of \epsilon',['N(0,',num2str(std(a(:),1),2),')'])
set(gca,'Title',text('String','Residual of L-R Deconvolution'))
print('hist_lucy.eps','-depsc2')