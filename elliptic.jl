using PlotlyJS
using Plots
using Unzip

function ϕ1(y)
    return cos(y)
end


function ϕ2(y)
    return 0
end


function ϕ3(x)
    return cos(x)
end


function ϕ4(x)
    return 0
end


function sol(data)
    return cos(data[1]) * cos(data[2])
end

function interpol(U, y, N, M)
    for i in 2:N-1
        for j in 2:M-1
            U[i, j] = (U[1, j] + (y[i] - y[1]) * 
                       (U[end, j] - U[1, j]) / (y[end] - y[1]))
        end
    end
end

function liebmann(x, y, hx, hy, ly, ε)
    N = length(x)
    M = length(y)

    U = zeros(length(x), length(y))
    # lets find ∂Ω
    for i in 1:M
        U[1, i] = ϕ1(y[i])
        U[end, i] = ϕ2(y[i])
    end
    for j in 1:N
        U[j, 1] = ϕ3(x[j])
        U[j, end] = ϕ4(x[j])
    end

    interpol(U, y, N, M)

    count = 0
    max_diff = ε + 1
    while max_diff > ε
        max_diff = 0
        for i in 2:N-1
            for j in 2:M-1
                new_val = (
                           (U[i + 1, j] * (hy^2) +
                             U[i - 1, j] * (hy^2) +
                             U[i, j + 1] * (hx^2) +
                             U[i, j - 1] * (hx^2)) / (2*hy^2 + 2*hx^2 - 2*hx^2 * hy^2)
                          )
                diff = abs(new_val - U[i, j])
                max_diff = max(max_diff, diff)
                U[i, j] = new_val
            end
        end
        count += 1
    end
    return U, count
end


function relax(x, y, hx, hy, ly, ω, ε)
    N = length(x)
    M = length(y)

    U = zeros(length(x), length(y))
    # lets find ∂Ω
    for i in 1:M
        U[1, i] = ϕ1(y[i])
        U[end, i] = ϕ2(y[i])
    end
    for j in 1:N
        U[j, 1] = ϕ3(x[j])
        U[j, end] = ϕ4(x[j])
    end

    interpol(U, y, N, M)

    count = 0
    max_diff = ε + 1
    while max_diff > ε
        max_diff = 0
        U_old = copy(U)
        for i in 2:N-1
            for j in 2:M-1
                U[i, j] = (
                           (1 - ω) * U_old[i, j] + 
                       ω * (U[i + 1, j] * (hy^2) +
                             U[i - 1, j] * (hy^2) +
                             U[i, j + 1] * (hx^2) +
                             U[i, j - 1] * (hx^2)) / 
                           (2*hy^2 + 2*hx^2 - 2*hx^2 * hy^2)
                          )
                diff = abs(U[i, j] - U_old[i, j])
                max_diff = max(max_diff, diff)
            end
        end
        count += 1
    end
    return U, count
end

function get_errors(U, U2, U3, N, M)
    err = 0
    err2 = 0
    for i in 1:N
        for j in 1:M
            err += (U[i, j] - U2[i, j])^2
            err2 += (U[i, j] - U3[i, j])^2
        end
    end
    return err/((N+1)*(M+1)), err2/((N+1)*(M+1))
end

function nm(A, B)
    return maximum(abs.(A - B))
end

function step_error(lx, ly, N, M, ω, ε)
    hs = []
    err_h = []
    # variation of h
    for i in range(10, 50, step=10)
        hx = lx / (i-1)*2
        hy = ly / (i-1)*2
        x = range(0, lx, step=hx)
        y = range(0, ly, step=hy)
        mesh = collect(Iterators.product(x, y))
        U = sol.(mesh)
        U2, count = liebmann(x, y, hx, hy, ly, ε)
        U3, count_relax = relax(x, y, hx, hy, ly, ω, ε)
        push!(err_h, (nm(U, U2), nm(U, U3)))
        push!(hs, hx)
        @info "count", count
        @info "count_relax", count_relax
    end
    return err_h, hs
end

lx = pi/2
ly = pi/2
N = 20 
M = 20

hx = lx / (N-1)
hy = ly / (M-1)

x = range(0, lx, step=hx)
y = range(0, ly, step=hy)

ε = 1e-10
ω = 1.2

mesh = collect(Iterators.product(x, y))
U = sol.(mesh)

U2, count = liebmann(x, y, hx, hy, ly, ε)

U3, count_relax = relax(x, y, hx, hy, ly, ω, ε)

srf2 = Plots.surface(x, y, U2, xlabel="x", ylabel="y")

# Plots.savefig(srf2, "srf2.png")

srf3 = Plots.surface(x, y, U3, xlabel="x", ylabel="y", c=:jet)

# Plots.savefig(srf3, "srf3.png")

er1, er2 = get_errors(U, U2, U3, N, M)

err_from_h, hs = step_error(lx, ly, N, M, ω, ε)

E1 = Plots.plot(hs, [getfield.(err_from_h, 1)], labels=["liebman"], 
                title="график погрешности от шага")
E2 = Plots.plot(hs, [getfield.(err_from_h, 2)], labels=["relax"],
                title="график погрешности от шага")

# Plots.savefig(E1, "err_liebman.png")
# Plots.savefig(E2, "err_relax.png")
