using PyPlot

folder = "ldc7"
get_file = var -> string("temp/sol/", folder, "/", var, ".csv")
# get_file = var -> string("temp/sol/", var, ".csv")
x = readcsv(get_file("x"))
y = readcsv(get_file("y"))
ux = readcsv(get_file("ux"))
uy = readcsv(get_file("uy"))
# kp = readcsv(get_file("kp"))

figure();
# pcolor(x, y, kp, cmap = "Greys");
quiver(x, y, ux, uy, ux, scale = 40);
axis([-1,1,-1,1]);

# figure();
# pcolor(x, y, kp, cmap = "Greys");
# axis([-1,1,-1,1]);

figure();
streamplot(x, y, ux, uy, density = 2, color = ux);
axis([-1,1,-1,1])
