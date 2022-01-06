function [x, info] = shermanMorrisonIteration(g,A,Z0)
%
% function [s, info] = shermanMorrisonIteration(g,A,Z0)
%
% Authors:
% (c) Julianne Chung (e-mail: jmchung@vt.edu)
% and Matthias Chung (e-mail: mcchung@vt.edu)
% and Joseph Slagel (email:slagelj@vt.edu) in December 2014
%
% Description:
% Given a least squares problem
% min x j j Ax - b j j 2ˆ2 + lambdaˆ2jjL x j j 2ˆ2
% with a invertible matrix Z 0=1/lambdaˆ2*inv(L'*L), this methods solves
% the normal equations
% (A'A + inv(Z0))x = g
% Where g=A'b. This method computes the vector
% x = inv((A'A + inv(Z0)))g
% by using the Sherman-Morrison formula and avoiding to build any matrix.
%
% Input arguments:
% g - A'b
% A - the A matrix in system Ax-b
% Z0 - Z0=1/lambdaˆ2*inv(L'L) where L is regularization matrix
%
% Output arguments:
% x - Solution to regularized least squares problem
% info - #
% .X - step approximations
[l,n] = size(A); % get size of A
x = Z0*g; % initialize solution vector
Z = Z0*A'; % initialize Z_{i,j} vectors
if nargout > 1, info.X = zeros(n,l); end
for i = 1:l % loop over all rows of A
    ff = 1+A(i,:)*Z(:,i);
    Z(:,i+1:l) = Z(:,i+1:l) - bsxfun(@times,A(i,:)*Z(:,i+1:l)/(1+A(i,:)*Z(:,i)),Z(:,i));
    %Update Z_{i,j} vectors
    x = x - bsxfun(@times,A(i,:)*x/(1+A(i,:)*Z(:,i)),Z(:,i));
    %Update x via Sherman Morisson formula
    if nargout > 1, info.X(:,i) = x; end
end
if nargout > 1, info.Z = Z; end