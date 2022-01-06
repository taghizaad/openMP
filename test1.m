clear
clc
row = 21;
col = 21;
A = normrnd(0,1,[row,col]);
lambda = 2;
L = eye(col);
Z0 = 1/lambda^2 * inv(L'*L);
U = A';
V = A';
x = ones(col,1);
b = A*x;
g = A'*b;

DesiredSD = 0.01 ;                               % the desired standard deviation
noise = random('ncx2',10*rand,10*rand,[col 1]) ;  % some random noise
noise = DesiredSD * (noise ./ std(noise)) ;        % scale the standard deviation
b_noisy = b + noise;

% x_ls_noisy = inv(A'*A)*A'*b_noisy;

x_ls = inv(A'*A)*g;

[x_tik, info] = shermanMorrisonIteration(g,A,Z0);

norm(x_tik - x_ls)
