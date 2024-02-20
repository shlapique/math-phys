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
            a = zeros(N)
            b = zeros(N)
            c = zeros(N)
            d = zeros(N)

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


N_x = 20
lx = 0
rx = 1

N_y = 20
ly = 0
ry = 1

h1 = (rx - lx) / N_x # x step
h2 = (rx - lx) / N_x # y step

K = 50
T = 1
τ = T / K # t step

α = 1

x = range(lx, rx, step=h1)
println("x:", x)

y = range(ly, ry, step=h2)
println("y:", y)

t = range(0, T, step=τ)
println("t:", t)


mesh = collect(Iterators.product(x, y, t))
U = sol.(mesh)
U2 = mvd(x, y, t, h1, h2, τ, α)


plt = Plots.plot(y, [U[10, :, 49] U2[10, :, 49]], labels=["точное решение" "mvd"])
