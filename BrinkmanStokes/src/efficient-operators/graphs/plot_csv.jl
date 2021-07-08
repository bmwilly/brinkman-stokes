using PyPlot

x = readcsv("/Users/bwilliams/Documents/brinkman-stokes/julia-parallel/efficient-operators/temp/sol/x.csv")
y = readcsv("/Users/bwilliams/Documents/brinkman-stokes/julia-parallel/efficient-operators/temp/sol/y.csv")
ux = readcsv("/Users/bwilliams/Documents/brinkman-stokes/julia-parallel/efficient-operators/temp/sol/ux.csv")
uy = readcsv("/Users/bwilliams/Documents/brinkman-stokes/julia-parallel/efficient-operators/temp/sol/uy.csv")

figure();
streamplot(x, y, ux, uy, density=4, color=ux);
axis([-1,1,-1,1]);

figure();
# pcolor(x, y, kp, cmap = "Greys");
quiver(x, y, ux, uy, ux, scale=2);
axis([-1,1,-1,1]);
