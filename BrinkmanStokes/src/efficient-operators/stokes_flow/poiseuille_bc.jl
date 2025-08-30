###POISEUILLE_BC specifies streamfunction associated with Poiseuille flow
# input
#   xbd       x boundary coordinate vector
#   ybd       y boundary coordinate vector
function poiseuille_bc(xbd, ybd)
	bc = ybd .* (1 - ybd .* ybd / 3) + (2 / 3)
	return bc
end
