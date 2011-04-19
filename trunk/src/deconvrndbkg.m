function [O_est,O_src]=deconvrndbkg(IMG,PSF,RES,NUMIT)
% DECONVRNDBKG deconvolution routine regularized by random background
% constraints (semi-monte-carlo approach).
% O_opt - Optimal estimate.
% O_est - Estimated Object.
% R_est - Estimation (Deconvolution) deviation.
% IMG - Observed Image.
% PSF - A priori PSF.
% RES - A priori Residual map used to extract background knowledge.
% NUMIT - Number of iterations for internal Bayesian deconvolution method.
NUMBG=numel(RES);
RES=RES(:);
O_est=zeros(numel(IMG),NUMBG);
O_src=O_est;
NPIXELS=numel(IMG);
tic
parfor n=1:NUMBG
    O_est(:,n)=reshape(deconvmap(IMG,PSF,RES(n),NUMIT),NPIXELS,[]);
    O_src(:,n)=O_est(:,n)-RES(n);
end
toc
save(['mc_',num2str(numel(RES)),'x',num2str(NUMIT)])
return