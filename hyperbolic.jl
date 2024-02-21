using PlotlyJS
using Plots
using Unzip

function ϕ0(t)
    return exp(-t) * cos(2*t)
end


function ϕl(t)
    return 0
end


function ψ1(x)
    return exp(-x) * cos(x)
end


function ψ2(x)
    return -exp(-x) * cos(x)
end


function sol(data)
    return exp(-data[2] - data[1]) * cos(data[1]) * cos(2*data[2])
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
    for i in 1:length(x)
        U[i, 1] = ψ1(x[i])
    end
    U[1, 2] = ϕ0(t[1])
    U[end, 2] = ϕl(t[end])
    # 1st order
    # for i in 1:length(x)
    #     U[i, 2] = ψ1(x[i]) + ψ2(x[i])*τ
    # end

    # 2nd order
    for i in 2:length(x)
        U[i, 2] = (
                   ψ1(x[i]) + ψ2(x[i])*τ 
                    + (2*exp(-x[i])*cos(x[i]) 
                    + 2*sin(x[i])/exp(x[i]) 
                    - (cos(x[i]) + sin(x[i]))/exp(x[i]) 
                    - 3*U[i, 1])*τ^2/2
                   )
    end

    for j in 3:length(t)
        for i in 2:length(x)-1
            U[i, j] = (1 / (h^2 * (2*τ + 1)) * (2*h^2 * (τ + 1) * U[i, j-1] 
                + (-3*U[i, j-1] * h^2 + (U[i+1, j-1] - U[i-1, j-1]) 
                * h - 2*U[i, j-1] 
                + U[i+1, j-1] 
                + U[i-1, j-1])*τ^2 
                - U[i, j-2]*h^2))
        end
        U[1, j] = ϕ0(t[j])
        U[end, j] = ϕl(t[j])
    end
    return U
end

function implicit(x, t, h, τ)
    U = zeros(length(x), length(t))
    σ = τ^2 / h^2
    N = length(x)
    K = length(t)

    # inital cond
    for i in 1:N
        U[i, 1] = ψ1(x[i])
    end
    U[1, 2] = ϕ0(t[1])
    U[end, 2] = ϕl(t[end])

    # inital cond [speed]
    # 1st order
    # for i in 1:N
    #     U[i, 2] = ψ1(x[i]) + ψ2(x[i])*τ
    # end

    # 2nd order
    for i in 2:N
        U[i, 2] = (
                   ψ1(x[i]) + ψ2(x[i])*τ 
                    + (2*exp(-x[i])*cos(x[i]) 
                    + 2*sin(x[i])/exp(x[i]) 
                    - (cos(x[i]) + sin(x[i]))/exp(x[i]) 
                    - 3*U[i, 1])*τ^2/2
                   )
    end

	for j in 3:K
        a = zeros(N-2)
        b = zeros(N-2)
        c = zeros(N-2)
        d = zeros(N-2)
        
        b[1] = -h^2 * (1 + 2*τ) - τ^2 * (2 + 3*h^2)
        c[1] = τ^2 * (1 + h)
        d[1] = h^2 * U[2, j-2] - (2 + 2*τ) * h^2 * U[2, j-1] - τ^2 * (1 - h) * ϕ0(t[j])
        a[end] = τ^2 * (1 - h)
        b[end] = -h^2 * (1 + 2*τ) - τ^2 * (2 + 3*h^2)
        d[end] = h^2 * U[N-1, j-2] - (2 + 2*τ) * h^2 * U[N-1, j-1] - τ^2 * (1 + h) * ϕl(t[j])
        for i in 3:N-2
            a[i-1] = τ^2 * (1 - h)
            b[i-1] = -h^2 * (1 + 2*τ) - τ^2 * (2 + 3*h^2)
            c[i-1] = τ^2 * (1 + h)
            d[i-1] = h^2 * U[i, j-2] - (2 + 2*τ) * h^2 * U[i, j-1]
        end
        res = tma(a, b, c, d)
        push!(res, ϕl(t[j]))
        insert!(res, 1, ϕ0(t[j]))
        U[:, j] = res
	end
	return U
end

function get_errors(U, U2, U3, N, K)
    err_explicit = 0
    err_implicit = 0
    for i in 1:N
        for j in 1:K
            err_explicit += (U[i, j] - U2[i, j])^2
            err_implicit += (U[i, j] - U3[i, j])^2
        end
    end
    return err_explicit/((N+1)*(K+1)), err_implicit/((N+1)*(K+1))
end

function nm(A, B)
    return maximum(abs.(A - B))
end

function step_error(l, T, N, K)
    # variation of τ
    err_τ = []
    taus = []
    hs = []
    err_h = []
    h_fixed = (l - 0) / N
    τ_fixed = T / K
    for i in range(100, 800, step=10)
        h = h_fixed
        x = range(0, l, step=h)
        τ = T / i
        t = range(0, T, step=τ)
        mesh = collect(Iterators.product(x, t))
        U = sol.(mesh)
        U2 = explicit(x, t, h, τ)
        U3 = implicit(x, t, h, τ)
        push!(err_τ, (nm(U, U2), nm(U, U3)))
        push!(taus, τ)
    end
    # variation of h
    for i in range(100, 800, step=10)
        h =  (l - 0) / i 
        x = range(0, l, step=h)
        τ = τ_fixed
        t = range(0, T, step=τ)
        mesh = collect(Iterators.product(x, t))
        U = sol.(mesh)
        U2 = explicit(x, t, h, τ)
        U3 = implicit(x, t, h, τ)
        push!(err_h, (nm(U, U2), nm(U, U3)))
        push!(hs, h)
    end
    return err_τ, taus, err_h, hs
end

T = 1
l = pi/2

a = 1

N = 1000   # segments for x
K = 1000  # segments for t

τ = T / K       # t step
println("τ=", τ)
h = (l - 0) / N # x step
println("h=", h)

# for explicit solution
σ = a^2 * τ^2 / h^2
println("σ=", σ)

x = range(0, l, step=h)
println("x:", x)
t = range(0, T, step=τ)
println("t:", t)

mesh = collect(Iterators.product(x, t))
U = sol.(mesh)

if σ <= 1
    println("Условие Куррента выполнено:", σ, "<= 1\n")
    U2 = explicit(x, t, h, τ)
    srf2 = Plots.surface(x, t, U2', c = :matter)
    plt2 = Plots.plot(x, [U[:, K] U2[:, K]], labels=["точное решение" "явная схема"])
else
    @warn "НЕ выполнено условие Куранта!"
end

srf2 = PlotlyJS.plot([PlotlyJS.surface(x=x, y=t, z=U2, name="explicit", colorscale="Jet"),
                      PlotlyJS.surface(x=x, y=t, z=U, name="solution")],
                     Layout(title="точное решение и явная схема"))
U3 = implicit(x, t, h, τ)

er1, er2 = get_errors(U, U2, U3, N, K)

srf3 = PlotlyJS.plot([PlotlyJS.surface(x=x, y=t, z=U3, name="implicit", colorscale="Blackbody"),
                      PlotlyJS.surface(x=x, y=t, z=U, name="solution")],
                     Layout(title="точное решение и неявная схема"))


err_from_τ, taus, err_from_h, hs = step_error(l, T, N, K)

E = Plots.plot(taus, [getfield.(err_from_τ, 1) getfield.(err_from_τ, 2)], labels=["явная" "неявная"], title="график погрешности от t")
E2 = Plots.plot(hs, [getfield.(err_from_h, 1) getfield.(err_from_h, 2)], labels=["явная" "неявная"], title="график погрешности от x")

# e1 = [nm(U[:, i], U2[:, i]) for i in 1:length(t)]
# e2 = [nm(U[:, i], U3[:, i]) for i in 1:length(t)]
# ep = Plots.plot(t, [e1 e2], labels=["явная" "неявная"], title="график погрешности от t")
