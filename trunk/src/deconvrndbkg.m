function [O_opt,O_est,R_est]=deconvrndbkg(I,P,R,NUMIT,NSIG)
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
NUMBG=length(R(:));
R=R(:);
O_est=zeros(size(I));
for n=1:NUMBG
    O_est=O_est+deconvmap(I,P,NUMIT,0,R(n));
end
O_est=O_est/NUMBG;
R_est=zeros(size(I));
for n=1:NUMBG
    R_est=(deconvmap(I,P,NUMIT,0,R(n))-O_est).^2+R_est;
end
R_est=sqrt(R_est/(NUMBG-1));
O_opt=double((O_est-mean(R))>=NSIG*R_est).*O_est+...
    double((O_est-mean(R))<NSIG*R_est)*mean(R);
return