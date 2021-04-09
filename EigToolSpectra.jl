using LinearAlgebra

#A is a matrix


funtion (A)

T = schur(A)

N = size(T)[1]

dx = -2:0.5:2
dy = -2:0.5:2

sigmin = zeros(Float64, length(dx), length(dy))
for (k, x) in enumerate(dx)
    for (j, y) in enumerate(dy)
    T1 = ((x+im*y)*I - T); T2 = transpose(T1);

    sigold = 0;qold = zeros(n,1);beta = 0;H = Float64[];

    q = randn(N,1)+im*randn(N,1); q = q/norm(q);

    for p in 1:maxit

    

        sigmin[j, k] = minimum((svd((x+im*y)*I-T)).S)
    end
end

end
