function result=mdmc(O,PSF,exposure,sig_obs,loopgain,maxnumloop,nummap,nummc)
%MDMC 1-D meta-deconvolution-monte-carlo experiment
% Simulating observation:
O=O(:);
PSF=PSF(:);
%I_obs=zeros(numel(O),nummc);
O_map=zeros(numel(O),nummc);
O_net=zeros(numel(O),nummc);
bg=zeros(1,nummc);
sig_bg=zeros(1,nummc);
numloop=zeros(1,nummc);
tic
parfor k=1:nummc
    I_obs=floor(imconv(normrnd(O*exposure,sqrt(O*exposure)),PSF)+...
        normrnd(0,sig_obs,size(O)));
    [~,rmap,sig_clean]=clean(I_obs,PSF,loopgain,maxnumloop);
    numloop(k)=sum(double(sig_clean>0));
    sig_bg(k)=sig_clean(numloop(k));
    if sig_clean(maxnumloop)~=0
        exit('max_num_loop for CLEAN algorithm is too low.')
    end
    bg(k)=mean(rmap);
    O_map(:,k)=deconvmap(I_obs,PSF,bg(k),nummap);
    O_net(:,k)=O_map(:,k)-bg(k);
end
toc
result=struct('bkg',bg,...
    'map',O_map,...
    'net',O_net,...
    'bkg_std',sig_bg,...
    'clean_num_loops',numloop);
save(['mdmc_exp=',num2str(numel(exposure)),...
    '_obsnoise=',num2str(sig_obs),...
    '_mcloop=',num2str(nummc),...
    '_map=',num2str(nummap)],'-v7.3')
return