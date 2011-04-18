function result=deconvmca(O_src,nbins)
%DECONVMCA DECONVolution Monte Carlo Analysis routine.
% O_src - Estimates of net objects (with the backgrounds subtracted).
dim=size(O_src);
npts=dim(length(dim));
npixels=numel(O_src)/npts;
O_src=reshape(O_src,npixels,[]);
o_mpll=zeros(npixels,1); % lower limit of most-probable interval.
o_mpul=zeros(npixels,1); % upper limit of most-probable interval.
%Prob_mp=zeros(npixels,1); % probability of most-probable interval.
prob_pos=zeros(npixels,1); % probability of positive interval.
%prob_ints=zeros(npixels,nbins+1); % probability of intervals.
parfor k=1:npixels
%    [o_mpll(k),o_mpul(k),prob_ints(k,:)]=mostprob(O_src(k,:),nbins);
    [o_mpll(k),o_mpul(k),~]=mostprob(O_src(k,:),nbins);
    prob_pos(k)=sum(O_src(k,:)>0)/npts;
end
%Prob_mp=Prob_mp/npts;
if length(dim)>2
    o_mpll=reshape(o_mpll,dim(1),[]);
    o_mpul=reshape(o_mpul,dim(1),[]);
    prob_pos=reshape(prob_pos,dim(1),[]);
%    prob_ints=reshape(prob_ints,dim(1),dim(2),[]);
end
result=struct('o_mpll',o_mpll,'o_mpul',o_mpul,'prob_pos',prob_pos);
return