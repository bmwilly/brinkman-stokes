
order = 2; msize = 4; dim = 2; 
t = hos_homg(order, msize, dim);

orders = [4];
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
