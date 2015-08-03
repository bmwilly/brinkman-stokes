module Tensor


function IAX(A::Array, x::Array)
    N = size(A, 1)
    y = A * reshape(x, N, N)
    y = y[:]
end

function AIX(A::Array, x::Array)
    N = size(A, 1)
    y = A * reshape(x, N, N)'
    y = y'
    y = y[:]
end

function IIAX(A::Array, x::Array)
    N = size(A, 1)
    y = A * reshape(x, N, N*N)
    y = y[:]
end
function IAIX(A::Array, x::Array)
    N = size(A, 1)
    q = reshape(x, N, N, N)
    y = zeros(N,N,N)
    for i=1:N
        y[i,:,:] = reshape(A * squeeze(q[i,:,:],1),1,N,N)
    end
    y = reshape(y,size(y)[1].*size(y)[2].*size(y)[3],1)
end
function AIIX(A::Array, x::Array)
    N = size (A, 1)
    y = reshape(x, N*N, N) * A'
    y = y[:]
end
function grad(refel, u)
    du = zeros(length(u), refel.dim);
    if (refel.dim == 2)
        du[:,1] = IAX(refel.Dr, u);
        du[:,2] = AIX(refel.Dr, u);
    else
        du[:,1] = IIAX(refel.Dr, u);
        du[:,2] = IAIX(refel.Dr, u);
        du[:,3] = AIIX(refel.Dr, u);
    end
    return du
end    
function grad2(A::Array, x::Array)
    IAX(A, x), AIX(A, x)
end
function grad3(A::Array, x::Array)
    IIAX(A, x), IAIX(A, x), AIIX(A, x)
end

end

