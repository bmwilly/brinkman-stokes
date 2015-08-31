% add lightspeed to path
addpath(genpath('/Users/bwilliams/Documents/masters-thesis/speed-test/matlab/lightspeed'))
addpath(genpath('/home/bmw313/Documents/masters-thesis/speed-test/matlab/lightspeed'))

% add files to path
addpath(genpath('/Users/bwilliams/Documents/masters-thesis/speed-test/matlab/global-matrix'))
% addpath(genpath('/home/bmw313/Documents/masters-thesis/speed-test/matlab/global-matrix'))

solve_square_stokes

tic;
for cnt = 1:100
  u = rand(nu + np, 1);
  w = K * u;
end
etoc = toc;

fprintf('vectors generated in %8.3e seconds\n',etoc)

numFlops = flops_mul(K, u);
fprintf('# of flops: %8.3e \n', numFlops)
