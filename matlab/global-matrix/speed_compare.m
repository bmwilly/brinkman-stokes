% orders = [2,3,4,5,6,7,8,9,10];
% orders = [11,12,13,14,15,16,17,18,19,20];
orders = [24];
msizes = [4];

otimes = [];
% dim = input('Dimension: ');
dim = 2;
if dim == 2
  msize = 5;
elseif dim == 3
  msize = 3;
end
for order = orders
  t = hos_homg(order, msize, dim);
  otimes = [otimes t];
end
