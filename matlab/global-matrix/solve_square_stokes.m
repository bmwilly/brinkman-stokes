% solve_square_stokes
% This sets up and solves the Stokes problem in a square domain,
% with a Q2-P1 discretization.
% It is known as driven-cavity flow, a model of the
% flow in a square cavity with the lid moving from left to right.
% It has a zero Dirichlet boundary condition.
% We use a regularized cavity model (as opposed to leaky or watertight).

% add files to path
addpath(genpath('/Users/bwilliams/Documents/masters-thesis/speed-test/matlab/global-matrix'))
addpath(genpath('/home/bmw313/Documents/masters-thesis/speed-test/matlab/global-matrix'))

% sets up problem
square_stokes

% solves problem
solve_stokes

% write results to file
% csvwrite('/home/bmw313/Documents/masters-thesis/speed-test/matlab/results/rslts.csv', )
