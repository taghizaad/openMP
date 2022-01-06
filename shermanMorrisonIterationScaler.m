function [x, info, U] = shermanMorrisonIterationScaler(g,A,lambda)
%
% function [s, info] = shermanMorrisonIteration(g,A,Z0)
%
% Authors:
% (c) Julianne Chung (e-mail: jmchung@vt.edu)
% and Matthias Chung (e-mail: mcchung@vt.edu)
% and Joseph Slagel (email:slagelj@vt.edu) in December 2014
%
% Description:
% Given a least squares problem min_x ||Ax - b||_2ˆ2 + lambdaˆ2||x||_2ˆ2
% This methods solves the normal equations (A'A + 1/lambdaˆ2*I)x = g
% Where I is the identity matrix and g=A'*b This method computes the vector x = inv(A'A + 1/lambdaˆ2*I)g
% by using the Sherman-Morrison formula and avoiding to build any matrix.
%
% Input arguments:
% g - A'b
% A - the A matrix in system Ax-b
% lambda - regularization paramater
% Output arguments:
% x - Solution to regularized least squares problem
% info - #
% .X - step approximations
lambda=1/lambda^2;
[l,n] = size(A); % get size of A
x = lambda*g; % initialize solution vector
U=lambda*A'; % initialize Z_{i,j} vectors
if nargout > 1, info.X = zeros(n,l); end
for i = 1:l % loop over all rows of A
    U(:,i+1:l) = U(:,i+1:l) - bsxfun(@times,A(i,:)*U(:,i+1:l)/(1+A(i,:)*U(:,i)),U(:,i));
    %Update Z_{i,j} vectors
    x = x - (A(i,:)*x)/(1+A(i,:)*U(:,i))*U(:,i);
    %Update x via Sherman Morisson formula
    if nargout > 1, info.X(:,i) = x; end
end
