# using PlotlyJS
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
function ϕ1(t)
    return -exp(-0.5 * t)
end

# U(x, t) = exp(-0,5t)sin(x)
# function sol(x, t)
#     return exp(-0.5 * t) * sin.(x)
# end
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
            U[i, j] = U[i, j-1] + (τ / h^2) * (U[i-1, j-1] - 2*U[i, j-1] + U[i+1, j-1]) + τ*f(x[i], t[j])
        end
        U[1, j] = U[2, j] - h * ϕ0(t[j])
        U[end, j] = h * ϕ1(t[j]) + U[end-1, j]
    end
    return U
end

function implicit_crank(x, t, h, τ, θ)
    if θ == 1
        println("Running implicit scheme...")
    else
        println("Running Crank–Nicolson method...")
    end
    U = zeros(length(x), length(t))
    σ = τ / h^2
    N = length(x)
    K = length(t)

    # inital cond
    for i in 1:N
        U[i, 1] = ψ(x[i])
    end
	for j in 2:K
        # println("ITER:", j)
        # display(U)
        a = zeros(N)
        b = zeros(N)
        c = zeros(N)
        d = zeros(N)
        
        # 2p approximation with 1st order
        b[1] = -1
        c[1] = 1
        d[1] = ϕ0(t[j])*h
        a[end] = -1
        b[end] = 1
        d[end] = ϕ1(t[j])*h
        # 3p approximation with 2nd order
        # b[1] = 2/h + h/τ
        # c[1] = -2/h
        # d[1] = h*U[1, j-1]/τ - ϕ0(t[j])*2
        # a[end] = -2/h
        # b[end] = 2/h + h/τ
        # d[end] = h*U[end, j-1] / τ+ϕ1(t[j])*2

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

function get_errors(U, U2, U3, N, K)
    err_explicit = 0
    err_implicit_crank = 0
    for i in 1:N
        for j in 1:K
            err_explicit += (U[i, j] - U2[i, j])^2
            err_implicit_crank += (U[i, j] - U3[i, j])^2
        end
    end
    return err_explicit/((N+1)*(K+1)), err_implicit_crank/((N+1)*(K+1))
end

T = 3
l = pi

a = 1


N = 50   # segments for x
K = 5000 # segments for t

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
t = range(0, T, step=τ)
print("t:")
println(t)

# println("Enter t:")
# t_cur = readline()
# t_cur = parse(Float64, t_cur)

t_cur = K*τ
mesh = collect(Iterators.product(x, t))
U = sol.(mesh)
plt = Plots.plot!(x, U[:, K])

if a * τ / h^2 <= 0.5
    println("Условие Куррента выполнено:", a * τ / h^2, "<= 0.5\n")
    U2 = explicit(x, t, h, τ)
    U3 = implicit_crank(x, t, h, τ, 1)

    srf2 = Plots.surface(x, t, U2')
    plt2 = Plots.plot(x, U2[:, K])
    plt3 = Plots.plot!(x, U3[:, K])
end
srf = PlotlyJS.plot(PlotlyJS.surface(x=x, y=t, z=U))

# get errors:
er1, er2 = get_errors(U, U2, U3, N, K)
println(er1, er2)
