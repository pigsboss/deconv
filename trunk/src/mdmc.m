function result=mdmc(O,PSF,exposure,sig_obs,loopgain,maxnumloop,nummap,nummc)
%MDMC 1-D meta-deconvolution-monte-carlo experiment
% Simulating observation:
O=O(:);
PSF=PSF(:);
%I_obs=zeros(numel(O),nummc);
O_map=zeros(numel(O),nummc);
bg=zeros(1,nummc);
parfor k=1:nummc
    I_obs=imconv(normrnd(O*exposure,sqrt(O*exposure)),PSF)+...
        normrnd(0,sig_obs,size(O));
    [~,rmap,sig_clean]=clean(I_obs,PSF,loopgain,maxnumloop);
    if sig_clean(maxnumloop)~=0
        exit('max_num_loop for CLEAN algorithm is too low.')
    end
    bg(k)=mean(rmap);
    O_map(:,k)=deconvmap(I_obs,PSF,bg(k),nummap);
end
result=struct('backgrounds',bg,...
    'map',O_map);
save(['mdmc_exp=',num2str(numel(exposure)),...
    '_obsnoise=',num2str(sig_obs),...
    '_mcloop=',num2str(nummc),...
    '_map=',num2str(nummap)],'-v7.3')
return