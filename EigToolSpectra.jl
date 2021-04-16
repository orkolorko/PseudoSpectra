using LinearAlgebra

function pseudospectr_Treph(A::Matrix, maxit = 10)
    #A is a matrix
    T = schur(A).T

    N = size(T)[1]

    dx = -2:0.5:2
    dy = -2:0.5:2

    sigmin = zeros(Float64, length(dx), length(dy))
    for (k, x) in enumerate(dx)
        for (j, y) in enumerate(dy)
            T1 = ((x+im*y)*I - T); T2 = transpose(T1);

            sigold = 0
            sig = 0
            qold = zeros(N,1);beta = 0
            H = Matrix{Float64}(undef, maxit+1, maxit+1)

            q = randn(N,1)+im*randn(N,1)
            q = q/norm(q)

            for p in 1:maxit
                v = T1\(T2\q) - beta*qold;

                alpha = real(q'*v)
                @info alpha[1]
                @info q
                v = v - alpha[1]*q

                beta = norm(v)
                qold = q
                q = v/beta

                H[p+1, p] = beta
                H[p, p+1] = beta
                H[p, p] = alpha[1]
                @H
                sig = maximum(eigvals(H[1:p, 1:p]))
                @info sig
                if abs(sigold/sig - 1)<1e-3
                     break
                end
                sigold  = sig;
            end
            @info sig
            sigmin[k,j] = sqrt(sig);
        end
    end
    return sigmin
end
