%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SVD Filter for multi-echo data sets
% Algorithm from: Bydder and Du (10.1016/j.mri.2006.03.006)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = SVD_filter(input)
[m n p] = size(input);
H = reshape(input,m.*n,p);
[U S V] = svd(H,0);
S = diag(S);
f = 1 - (S(p)./S).^2;
S = diag(f.*S);
H = U*S*V';
output = reshape(H,m,n,p);
