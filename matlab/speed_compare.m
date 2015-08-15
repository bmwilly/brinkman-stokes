orders = [4];
msizes = [4];

otimes = [];
dim = input("Dimension: ");
if dim == 2
  msize = 5;
elseif dim == 3
  msize = 3;
end
for order in orders
  t = hos_homg(order, msize, dim);
  otimes = [otimes t];
end
