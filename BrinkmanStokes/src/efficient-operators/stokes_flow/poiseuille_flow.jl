###POISEUILLE_FLOW specifies Poiseuille flow boundary condition
# input
#   xbd       x coordinate vector
#   ybd       y coordinate vector
function poiseuille_flow(xbd, ybd)
	bcx = 1 .- ybd .* ybd
	bcy = 0 .* xbd
	(bcx, bcy)
end
