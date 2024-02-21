using PlotlyJS
using Plots
using Unzip

function ϕ1()
    return 0
end

function ϕ2(y, t)
    return y * cos(t)
end

function ϕ3()
    return 0
end

function ϕ4(x, t)
    return x * cos(t)
end

function ψ(x, y)
    return x * y
end

function f(x, y, t)
    return -x * y * sin(t)
end

function sol(data)
    return data[1] * data[2] * cos(data[3])
end

function tma(a, b, c, d)
    n = length(a)
    x = copy(d)
    c_prime = copy(c)
    c_prime[1] /= b[1]
    x[1] /= b[1]
    for i in 2:n
        scale = 1.0 / (b[i] - c_prime[i-1]*a[i])
        c_prime[i] *= scale
        x[i] = (x[i] - a[i] * x[i-1]) * scale
    end
    # Back
    for i in n-1:-1:1
        x[i] -= (c_prime[i] * x[i+1])
    end

    return x
end


# method of variable directions
function mvd(x, y, t, h1, h2, τ, α)
    N = length(x)
    M = length(y)
    K = length(t)
    U = zeros(N, M, K)
    for i in 1:N
        for j in 1:M
            U[i, j, 1] = ψ(x[i], y[j])
        end
    end

    for k in 2:K
        # layer
        u1 = zeros(N, M)
        t_half = t[k] - τ / 2
        uk_1 = U[:, :, k-1]
        
        for j in 2:M-1
            a = zeros(N)
            b = zeros(N)
            c = zeros(N)
            d = zeros(N)
            
            b[1] = 1
            c[1] = 0
            d[1] = ϕ1()

            a[end] = 0
            b[end] = 1
            d[end] = ϕ2(y[j], t_half)

            for i in 2:N-1
                a[i] = α 
                b[i] = -2 * (h1^2) / τ - 2 * α
                c[i] = α 
                d[i] = (-(h1^2) * (uk_1[i, j+1] - 2 * uk_1[i, j] + uk_1[i, j-1]) / (h2^2) -
                        2 * (h1^2) * uk_1[i, j] / τ - (h1^2) * f(x[i], y[j], t_half))
            end

            xu = tma(a, b, c, d)

            for i in 1:N
                u1[i, j] = xu[i]
                u1[i, 1] = ϕ3()
                u1[i, end] = ϕ4(x[i], t_half)
            end
        end

        u2 = zeros(N, M)

        for i in 2:N-1
            a = zeros(M)
            b = zeros(M)
            c = zeros(M)
            d = zeros(M)

            b[1] = 1
            c[1] = 0
            d[1] = ϕ3()

            a[end] = 0
            b[end] = 1
            d[end] = ϕ4(x[i], t[k])

            for j in 2:M-1
                a[j] = α
                b[j] = -2 * (h2^2) / τ - 2 * α
                c[j] = α
                d[j] = (-2 * (h2^2) * u1[i, j] / τ - 
                α * (h2^2) * (u1[i+1, j] - 2 * u1[i, j] + 
                              u1[i-1, j]) / (h1^2) - (h2^2) * f(x[i], y[j], t_half))
            end

            xu = tma(a, b, c, d)

            for j in 1:M
                u2[i, j] = xu[j]
                u2[1, j] = ϕ1()
                u2[end, j] = ϕ2(y[j], t[k])
            end
        end

        for i in 1:N
            for j in 1:M
                U[i, j, k] = u2[i, j]
            end
        end
    end
    return U
end


# fractional step method
function fsm(x, y, t, h1, h2, τ, α)
    N = length(x)
    M = length(y)
    K = length(t)
    U = zeros(N, M, K)
    for i in 1:N
        for j in 1:M
            U[i, j, 1] = ψ(x[i], y[j])
        end
    end

    for k in 2:K
        # layer
        u1 = zeros(N, M)
        t_half = t[k] - τ / 2
        uk_1 = U[:, :, k-1]
        
        for j in 2:M-1
            a = zeros(N)
            b = zeros(N)
            c = zeros(N)
            d = zeros(N)
            
            b[1] = 1
            c[1] = 0
            d[1] = ϕ1()

            a[end] = 0
            b[end] = 1
            d[end] = ϕ2(y[j], t_half)

            for i in 2:N-1
                a[i] = α 
                b[i] = -(h1^2) / τ - 2 * α
                c[i] = α 
                d[i] = -(h1^2) * uk_1[i, j] / τ - (h1^2) * f(x[i], y[j], t_half) / 2
            end

            xu = tma(a, b, c, d)

            for i in 1:N
                u1[i, j] = xu[i]
                u1[i, 1] = ϕ3()
                u1[i, end] = ϕ4(x[i], t_half)
            end
        end

        u2 = zeros(N, M)

        for i in 2:N-1
            a = zeros(M)
            b = zeros(M)
            c = zeros(M)
            d = zeros(M)

            b[1] = 1
            c[1] = 0
            d[1] = ϕ3()

            a[end] = 0
            b[end] = 1
            d[end] = ϕ4(x[i], t[k])

            for j in 2:M-1
                a[j] = α
                b[j] = -(h2^2) / τ - 2 * α
                c[j] = α
                d[j] = -(h2^2) * u1[i, j] / τ - (h2^2) * f(x[i], y[j], t[k]) / 2
            end

            xu = tma(a, b, c, d)

            for j in 1:M
                u2[i, j] = xu[j]
                u2[1, j] = ϕ1()
                u2[end, j] = ϕ2(y[j], t[k])
            end
        end

        for i in 1:N
            for j in 1:M
                U[i, j, k] = u2[i, j]
            end
        end
    end
    return U
end


function get_errors(U, U2, U3, N, M)
    err_1 = 0
    err_2 = 0
    for i in 1:N
        for j in 1:M
            err_1 += (U[i, j] - U2[i, j])^2
            err_2 += (U[i, j] - U3[i, j])^2
        end
    end
    return err_1/((N+1)*(K+1)), err_2/((N+1)*(K+1))
end
function nm(A, B)
    return maximum(abs.(A - B))
end

function step_error(lx, rx, ly, ry, T, N, M, K, α)
    err_τ = []
    err_h1 = []
    err_h2 = []
    h1s = []
    h2s = []
    taus = []
    h1_fixed = (rx - lx) / N
    h2_fixed = (ry - lx) / M
    τ_fixed = T / K
    for i in range(10, 100, step=2)
        h1 = (rx - lx) / i
        h2 = h2_fixed
        τ = τ_fixed
        x = range(lx, rx, step=h1)
        y = range(ly, ry, step=h2)
        t = range(0, T, step=τ)
        mesh = collect(Iterators.product(x, y, t))
        U = sol.(mesh)
        U2 = mvd(x, y, t, h1, h2, τ, α)
        U3 = fsm(x, y, t, h1, h2, τ, α)
        push!(err_h1, (nm(U, U2), nm(U, U3)))
        push!(h1s, h1)
    end
    for i in range(10, 100, step=2)
        h1 = h1_fixed
        h2 = (rx - lx) / i
        τ = τ_fixed
        x = range(lx, rx, step=h1)
        y = range(ly, ry, step=h2)
        t = range(0, T, step=τ)
        mesh = collect(Iterators.product(x, y, t))
        U = sol.(mesh)
        U2 = mvd(x, y, t, h1, h2, τ, α)
        U3 = fsm(x, y, t, h1, h2, τ, α)
        push!(err_h2, (nm(U, U2), nm(U, U3)))
        push!(h2s, h2)
    end
    # variation of τ
    for i in range(10, 50, step=2)
        h1 = h1_fixed
        h2 = h2_fixed
        τ = T / i
        x = range(lx, rx, step=h1)
        y = range(ly, ry, step=h2)
        t = range(0, T, step=τ)
        mesh = collect(Iterators.product(x, y, t))
        U = sol.(mesh)
        U2 = mvd(x, y, t, h1, h2, τ, α)
        U3 = fsm(x, y, t, h1, h2, τ, α)
        push!(err_τ, (nm(U, U2), nm(U, U3)))
        push!(taus, τ)
    end
    return err_h1, h1s, err_h2, h2s, err_τ, taus
end

N = 20
lx = 0
rx = 1

M = 20
ly = 0
ry = 1

h1 = (rx - lx) / N # x step
h2 = (rx - lx) / M # y step

K = 50
T = 1
τ = T / K # t step

α = 1

t_test_index = 10
# x_test_index = 10
# y_test_index = 10


x = range(lx, rx, step=h1)
println("x:", x)

y = range(ly, ry, step=h2)
println("y:", y)

t = range(0, T, step=τ)
println("t:", t)


mesh = collect(Iterators.product(x, y, t))
U = sol.(mesh)
U2 = mvd(x, y, t, h1, h2, τ, α)
U3 = fsm(x, y, t, h1, h2, τ, α)

err_from_h1, h1s, err_from_h2, h2s, err_from_τ, taus = step_error(lx, rx, ly, ry, T, N, M, K, α)

E = Plots.plot(h1s, [getfield.(err_from_h1, 1) getfield.(err_from_h1, 2)], labels=["mvd" "fsm"], title="график погрешности от h1")
E2 = Plots.plot(h2s, [getfield.(err_from_h2, 1) getfield.(err_from_h2, 2)], labels=["mvd" "fsm"], title="график погрешности от h2")
E3 = Plots.plot(taus, [getfield.(err_from_τ, 1) getfield.(err_from_τ, 2)], labels=["mvd" "fsm"], title="график погрешности от t")

# plt = Plots.plot(y, [U[x_test_index, :, t_test_index] U2[x_test_index, :, t_test_index] U3[x_test_index, :, t_test_index]], 
# plt = Plots.plot(y, [U[x_test_index, :, t_test_index] U2[x_test_index, :, t_test_index] U3[x_test_index, :, t_test_index]], 
#                  labels=["точное решение" "mvd" "fsm"], 
#                  linestyle=[:solid :dash :dot], lw=3)

# er1, er2 = get_errors(U[:, :, t_test_index], U2[:, :, t_test_index], U3[:, :, t_test_index], N, M)


# e1 = [nm(U[:, i], U2[:, i]) for i in 1:length(t)]
# e2 = [nm(U[:, i], U3[:, i]) for i in 1:length(t)]
# ep = Plots.plot(t, [e1 e2], labels=["явная" "неявная"], title="график погрешности от t")

# e1_x = [nm(U[j, :], U2[j, :]) for j in 1:length(x)]
# e2_x = [nm(U[j, :], U3[j, :]) for j in 1:length(x)]
# ep_x = Plots.plot(x, [e1_x e2_x], labels=["явная" "неявная"], title="график погрешности от x")
