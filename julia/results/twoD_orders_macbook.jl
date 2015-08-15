using DataFrames
using Gadfly

df = readtable("results/twoD_orders_macbook.csv")
p = Gadfly.plot(x = df[:Order], y = df[:Time], color = df[:Method], Geom.line,
                Scale.x_discrete(),
                Guide.colorkey("Method"),
                Guide.xlabel("Order"), Guide.ylabel("Time (s)"))
draw(PNG("results/twoD_orders_macbook.png", 8inch, 4inch), p)
