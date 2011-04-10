function [ul,ll]=confint(x,X,level)
%CONFINT Confident interval estimate
% ul    - upper limit, scalar or row matrix
% ll    - lower limit, scalar or row matrix
% x     - point estimate, scalar or column vector
% X     - all samples, row vector or matrix
% level - confident level
X=sort(X,1);

return