using PyPlot

folder = "k4"
get_file = var -> string("temp/sol/", folder, "/", var, ".csv")
x = readcsv(get_file("x"))
y = readcsv(get_file("y"))
ux = readcsv(get_file("ux"))
uy = readcsv(get_file("uy"))

figure();
# pcolor(x, y, kp, cmap = "Greys");
quiver(x, y, ux, uy, ux, scale = 2);
axis([-1,1,-1,1]);
