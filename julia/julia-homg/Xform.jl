module Xform

function identity (Xin)
	Xout = copy(Xin)
end

function twoX (Xin)
	Xout = copy(Xin)
	Xout[:,1] = 2*Xout[:,1]
	return Xout
end

function twoY (Xin)
	Xout = copy(Xin)
	Xout[:,2] = 2*Xout[:,2];
	return Xout
end

function twoSix (Xin)
	# currently only 2D 
	Xout = copy(2*Xin)
	Xout[:,2] = 3*Xout[:,2];
	return Xout
end

function shell (Xin)
	d = size(Xin, 2);
	R2 = 1.0; # hard coded for now.
	R1 = 0.7; # hard coded for now.
	R2byR1 = R2 / R1;
	R1sqrbyR2 = R1 * R1 / R2;
	if (d == 2)
		x = zeros( size(Xin[:,1]) );
		y = tan ( Xin[:,1]  * pi/4 );
		R = R1sqrbyR2 * ( R2byR1.^(Xin[:,2] + 1) ) ;
	else
		x = tan ( Xin[:,1]  * pi/4 );
		y = tan ( Xin[:,2]  * pi/4 );
		R = R1sqrbyR2 * ( R2byR1.^(Xin[:,3] + 1) );
	end
	q = R ./ sqrt (x.*x + y.*y + 1);
	if (d == 3)
		Xout = zeros(size(q)[1],3) #define Xout to use below
		Xout[:,1] =  q.* y;
		Xout[:,2] = -q.* x;
		Xout[:,3] =  q;
	else
		Xout = zeros(size(q)[1],2) #define Xout to use below
		Xout[:,1] =   q.* y;
		Xout[:,2] =   q;
	end
	return Xout
end

end