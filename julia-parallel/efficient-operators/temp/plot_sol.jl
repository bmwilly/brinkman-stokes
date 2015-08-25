using PyPlot

x = readcsv("temp/sol/x.csv")
y = readcsv("temp/sol/y.csv")
ux = readcsv("temp/sol/ux.csv")
uy = readcsv("temp/sol/uy.csv")

figure();
# pcolor(x, y, kp, cmap = "Greys");
quiver(x, y, ux, uy, ux, scale = 2);
axis([-1,1,-1,1]);
