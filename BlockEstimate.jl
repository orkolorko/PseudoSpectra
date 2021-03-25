using LinearAlgebra
function startmatrix(n, t, T)
    X = 2*rand(n, t).-1
    norms = []
    for i in 1:t
        X[:,i]/=norm(X[:,i], 1)
    end
    return X
end

function canonicalbasis(i, n)
    e = zeros(n)
    e[i]+=1
    return e
end

function notparallel(S, Sold)
    t = size(S)[2]
    n = size(S)[1]
    iter = 0
    for i in 1:t
        while iter<n/t
            y = Sold'*S[:,i]
            iter +=1
            if norm(y, Inf)<n
                y = S[:,1:i-1]'*S[:,i]
                iter+=1
                if norm(y, Inf)<n
                    return S
                end
            end
        end
        S[:, i] = rand( (-1,1), n)
    end
    return S
end

function normestimate(A::Matrix{T}, t = 2, itmax = 3, rettype = :est) where {T}
    n = size(A)[1]
    X = startmatrix(n,t, T)
    ind_hist = Int64[]
    est_old = 0
    ind = zeros(n, 1)
    S = zeros(n, t)
    Sold = S
    ind_best = undef
    est = 0
    w = T[]
    for k in 1:itmax
        Y = A*X
        normY = [norm(Y[:,j], 1) for j in 1:t]
        #@info normY
        est, j = findmax(normY)
        #@info est_old
        #@info est, j
        if est>est_old || k==2
            ind_best = ind[j]
            w = Y[:, j]
        end
        if k>=2 && est<=est_old
            est = est_old
            @goto six
        end
        est_old = est
        S_old = S
        if k>itmax
            @goto six
        end
        S = sign.(Y)
        if t>1
            S = notparallel(S, Sold)
        end
        Z = A'*S
        h = [(i, norm(Z[i, :])) for i in 1:n]

        # please remark that I am using j as index variable again, it is not the same j
        esth, j = findmax(h)
        if k>=2 && ind_best==j
            @goto six
        end
        h =  sort!(h, lt = (x,y)->x[2]<y[2])
        ind = map(x->x[1], h)
        if t>1
            if issubset(ind, ind_hist)
                @goto six
            end
            ind = setdiff(ind, ind_hist) # check if setdiff preserves order, seems so
        end

        # please remark that I am using j as index variable again, it is not the same j
        for j in 1:t
            X[:, j] = canonicalbasis(ind[j],n)
        end
        append!(ind_hist,ind[1:t])
    end

    @label six
    v = canonicalbasis(ind_best, n)
    if rettype==:est
        return est
    end
    if rettype==:all
        return est, v, w
    end
end

function normestimate(A::Matrix{Complex{T}}, t = 2, itmax = 3, rettype = :est) where {T}
    n = size(A)[1]
    X = startmatrix(n,t, T)
    ind_hist = Int64[]
    est_old = 0
    ind = zeros(n, 1)
    S = zeros(n, t)
    Sold = S
    ind_best = undef
    est = 0
    w = T[]
    for k in 1:itmax
        Y = A*X
        normY = [norm(Y[:,j], 1) for j in 1:t]
        #@info normY
        est, j = findmax(normY)
        #@info est_old
        #@info est, j
        if est>est_old || k==2
            ind_best = ind[j]
            w = Y[:, j]
        end
        if k>=2 && est<=est_old
            est = est_old
            @goto six
        end
        est_old = est
        S_old = S
        if k>itmax
            @goto six
        end
        S = sign.(Y)
        Z = A'*S
        h = [(i, norm(Z[i, :])) for i in 1:n]

        # please remark that I am using j as index variable again, it is not the same j
        esth, j = findmax(h)
        if k>=2 && ind_best==j
            @goto six
        end
        h =  sort!(h, lt = (x,y)->x[2]<y[2])
        ind = map(x->x[1], h)
        if t>1
            if issubset(ind, ind_hist)
                @goto six
            end
            ind = setdiff(ind, ind_hist) # check if setdiff preserves order, seems so
        end

        # please remark that I am using j as index variable again, it is not the same j
        for j in 1:t
            X[:, j] = canonicalbasis(ind[j],n)
        end
        append!(ind_hist,ind[1:t])
    end

    @label six
    v = canonicalbasis(ind_best, n)
    if rettype==:est
        return est
    end
    if rettype==:all
        return est, v, w
    end
end
