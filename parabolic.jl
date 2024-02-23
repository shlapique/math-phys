using PlotlyJS
using Plots
using Unzip

# u(x, 0) = sin(x)
function ψ(x)
    return sin(x)
end

# u(0, t)
function ϕ0(t)
    return exp(-0.5 * t)
end

# u(pi, t)
function ϕl(t)
    return -exp(-0.5 * t)
end


function sol(data)
    return exp(-0.5 * data[2]) * sin(data[1])
end

function f(x, t)
    return 0.5 * exp(-0.5 * t) * sin(x)
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

function explicit(x, t, h, τ)
    U = zeros(length(x), length(t))

    # inital cond
    for i in 1:length(x)
        U[i, 1] = ψ(x[i])
    end

    for j in 2:length(t)
        for i in 2:length(x)-1
            U[i, j] = (U[i, j-1] + (τ / h^2) * (U[i-1, j-1] - 2*U[i, j-1] +
                                                U[i+1, j-1]) + τ*f(x[i], t[j]))
        end
        U[1, j] = U[2, j] - h * ϕ0(t[j])
        U[end, j] = h * ϕl(t[j]) + U[end-1, j]
    end
    return U
end

function implicit_crank(x, t, h, τ, θ)
    if θ == 1
        @info "Running implicit scheme..."
    else
        @info "Running Crank–Nicolson method..."
    end
    U = zeros(length(x), length(t))
    N = length(x)
    K = length(t)

    # inital cond
    for i in 1:N
        U[i, 1] = ψ(x[i])
    end
	for j in 2:K
        a = zeros(N)
        b = zeros(N)
        c = zeros(N)
        d = zeros(N)
        
        b[1] = -1
        c[1] = 1
        d[1] = ϕ0(t[j])*h
        a[end] = -1
        b[end] = 1
        d[end] = ϕl(t[j])*h

        for i in 2:N-1
            a[i] = θ*τ
            b[i] = -2*θ*τ - h^2
            c[i] = θ*τ
            d[i] = (
                    (θ-1)*τ*U[i-1, j-1] 
                    + (-h^2+2*τ-θ*τ*2)*U[i, j-1] 
                    + (θ-1)*τ*U[i+1, j-1] 
                    - θ*τ*h^2*f(x[i], t[j]) 
                    + (θ-1)*τ*h^2*f(x[i], t[j-1])
                   )
        end
        res = tma(a, b, c, d)
        U[:, j] = res 
	end
	return U
end

function get_errors(U, U2, U3, U3_crank, N, K)
    err_explicit = 0
    err_implicit = 0
    err_crank = 0
    for i in 1:N
        for j in 1:K
            err_explicit += (U[i, j] - U2[i, j])^2
            err_implicit += (U[i, j] - U3[i, j])^2
            err_crank += (U[i, j] - U3_crank[i, j])^2
        end
    end
    return err_explicit/((N+1)*(K+1)), err_implicit/((N+1)*(K+1)),
    err_crank/((N+1)*(K+1))
end

function nm(A, B)
    return maximum(abs.(A - B))
end

function step_error(l, T, N, K, σ)
    hs = []
    err_h = []
    # variation of h
    τ_fixed = T / K
    for i in range(10, 40, step=2)
        h =  (l - 0) / i 
        x = range(0, l, step=h)
        τ = τ_fixed
        t = range(0, T, step=τ)
        mesh = collect(Iterators.product(x, t))
        U = sol.(mesh)
        U2 = explicit(x, t, h, τ)
        U3 = implicit_crank(x, t, h, τ, 1)
        U3_crank = implicit_crank(x, t, h, τ, 0.5)
        push!(err_h, (nm(U, U2), nm(U, U3), nm(U, U3_crank)))
        push!(hs, h)
    end
    return err_h, hs
end

T = 1
l = pi

a = 1

N = 40   # segments for x
K = 500  # segments for t

θ = 1

τ = T / K       # t step
println("τ=", τ)
h = (l - 0) / N # x step
println("h=", h)

# for explicit solution
σ = a^2 * τ / h^2
println("σ=", σ)

x = range(0, l, step=h)
print("x:")
println(x)
@warn size(x)
t = range(0, T, step=τ)
print("t:")
println(t)

# println("Enter t:")
# t_cur = readline()
# t_cur = parse(Float64, t_cur)

t_cur = K*τ
mesh = collect(Iterators.product(x, t))
U = sol.(mesh)
plt = Plots.plot(x, U[:, K])

if σ <= 0.5
    println("Условие Куррента выполнено:", σ, "<= 0.5\n")
    U2 = explicit(x, t, h, τ)

    srf2 = Plots.surface(x, t, U2', xlabel="x", ylabel="time")
    # Plots.savefig(srf2, "srf2.png")
    plt2 = Plots.plot(x, [U[:, K] U2[:, K]], labels=["точное решение" "явная схема"])
end

U3 = implicit_crank(x, t, h, τ, 1)
srf3 = Plots.surface(x, t, U3', xlabel="x", ylabel="time", c=:buda)
plt3 = Plots.plot(x, [U[:, K] U3[:, K]], labels=["точное решение" "неявная схема"])
# Plots.savefig(srf3, "srf3.png")

U3_crank = implicit_crank(x, t, h, τ, 0.5)
srf_crank = Plots.surface(x, t, U3_crank', xlabel="x", ylabel="time", c=:jet)
# Plots.savefig(srf_crank, "srf_crank.png")
plt3_crank = Plots.plot(x, [U[:, K] U3_crank[:, K]], 
                        labels=["точное решение" "схема Кранка-Николсона"])

srf = PlotlyJS.plot(PlotlyJS.surface(x=x, y=t, z=U))

# get errors:
er1, er2, er3 = get_errors(U, U2, U3, U3_crank, N, K)
println("ERRORS:")
println(er1, "\n", er2, "\n", er3)

err_from_h, hs = step_error(l, T, N, K, σ)

E1 = Plots.plot(hs, [getfield.(err_from_h, 1)], 
               labels=["явная"], title="график погрешности от шага h")
E2 = Plots.plot(hs, [getfield.(err_from_h, 2)], 
               labels=["неявная"], title="график погрешности от шага h")
E3 = Plots.plot(hs, [getfield.(err_from_h, 3)], 
               labels=["Crank-Nicolson"], title="график погрешности от шага h")
# Plots.savefig(E1, "err_explicit.png")
# Plots.savefig(E2, "err_implicit.png")
# Plots.savefig(E3, "err_crank.png")
