using LinearAlgebra

#A is a matrix

T = schur(A)

N = size(T)[1]

dx = -2:0.5:2
dy = -2:0.5:2

sigmin = zeros(Float64, length(dx), length(dy))
for (k, x) in enumerate(dx)
    for (j, y) in enumerate(dy)
    T1 = ((x+im*y)*I - T); T2 = transpose(T1);

    sigold = 0;qold = zeros(N,1);beta = 0;H = Float64[];

    q = randn(N,1)+im*randn(N,1); q = q/norm(q);

    for p in 1:maxit
        v = T1\(T2\q) - beta*qold;

        alpha = real(Transpose(q)*v);v = v - alpha*q;

        beta = norm(v);qold = q;q = v/beta;

        H[p+1,p] = beta;H[p,p+1] = beta;H[p,p] = alpha;
        sig = max(eigvals(H[1:p,1:p]));
        if abs(sigold/sig - 1)<1e-3, break,end
        sigold  = sig;
    end

    sigmin[k,j] = sqrt(sig);

end

end
