function [O_avg,O_est,R_est]=deconvrndbkg(I,P,R,NUMIT)
% DECONVRNDBKG deconvolution routine regularized by random background
% constraints (semi-monte-carlo approach).
% O_opt - Optimal estimate.
% O_est - Estimated Object.
% R_est - Estimation (Deconvolution) deviation.
% I - Observed Image.
% P - A priori PSF.
% R - A priori Residual map used to extract background knowledge.
% NUMIT - Number of iterations for internal Bayesian deconvolution method.
% NSIG - N*sigma significance.
NUMBG=numel(R);
R=R(:);
O_est=zeros(numel(I),NUMBG);
for n=1:NUMBG
    O_est(:,n)=reshape(deconvmap(I,P,NUMIT,R(n)),numel(I),[]);
end
O_avg=mean(O_est,2);
R_est=std(size(I),1,2);
% O_opt=double((O_est-mean(R))>=NSIG*R_est).*O_est+double((O_est-mean(R))<NSIG*R_est)*mean(R);
return