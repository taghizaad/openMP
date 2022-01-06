function [x, info] = smi(b,U,V,Z0)

[t,k] = size(V); % get size of U
x = Z0*b; % initialize solution vector
Z = Z0*U; % initialize Z_{i,j} vectors
if nargout > 1, info.X = zeros(t,k); end
for i = 1:k % loop over all cols of V
    Z(:,i+1:k) = Z(:,i+1:k) - bsxfun(@times,V(:,i)'*Z(:,i+1:k)/(1+V(:,i)'*Z(:,i)),Z(:,i));
    %Update Z_{i,j} vectors
    x = x - bsxfun(@times,V(:,i)'*x/(1+V(:,i)'*Z(:,i)),Z(:,i));
    %Update x via Sherman Morisson formula
    if nargout > 1, info.X(:,i) = x; end
end
if nargout > 1, info.Z = Z; end