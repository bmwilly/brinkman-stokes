###ITERAPP apply matrix operator to vector and error gracefully
# Applies matrix operator AFUN to vector X.
# If ATYPE is "matrix", the AFUN is a matrix and the OP is applied directly.
# OP is either "mtimes" or "mldivide".
# ATYPE and AFCNSTR are used in case of error.
# ITERAPP is designed for use by iterative methods.
function iterapp(op, afun, atype, afcnstr, x, varargin)

	Af =  isa(A,Function) ? A : x->A*x
	# if isequal(atype, "matrix")
	if isa(atype, Function)
		if (nargin >= 6) & isequal(varargin{end}, "notransp")

		end
	end
end