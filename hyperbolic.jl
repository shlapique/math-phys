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
        # U[i, 2] = ψ1(x[i]) + ψ2(x[i])*τ + exp(-x[i])*sin(x[i])*τ^2
        U[i, 2] = ψ1(x[i]) + ψ2(x[i])*τ + (2*exp(-x[i])*cos(x[i]) + 2*sin(x[i])/exp(x[i]) - (cos(x[i]) + sin(x[i]))/exp(x[i]) - 3*U[i, 1])*τ^2/2
    end

    for j in 3:length(t)
        for i in 2:length(x)-1
            U[i, j] = 1 / (h^2 * (2*τ + 1)) * 
            (2*h^2 * (τ + 1) * U[i, j-1] + 
             (-3*U[i, j-1] * h^2 + (U[i+1, j-1] - U[i-1, j-1]) * h - 2*U[i, j-1] + U[i+1, j-1] + U[i-1, j-1])*τ^2 - U[i, j-2]*h^2)
        end
        U[1, j] = ϕ0(t[j])
        U[end, j] = ϕl(t[j])
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

# function get_errors(U, U2, U3, N, K)
#     err_explicit = 0
#     err_implicit_crank = 0
#     for i in 1:N
#         for j in 1:K
#             err_explicit += (U[i, j] - U2[i, j])^2
#             err_implicit_crank += (U[i, j] - U3[i, j])^2
#         end
#     end
#     return err_explicit/((N+1)*(K+1)), err_implicit_crank/((N+1)*(K+1))
# end
function get_errors(U, U2, N, K)
    err_explicit = 0
    for i in 1:N
        for j in 1:K
            err_explicit += (U[i, j] - U2[i, j])^2
        end
    end
    return err_explicit/((N+1)*(K+1))
end

function max_abs_error(A, B)
    return maximum(abs.(A - B))
end

T = 1
l = pi/2

a = 1

N = 40   # segments for x
K = 500  # segments for t

τ = T / K       # t step
println("τ=", τ)
h = (l - 0) / N # x step
println("h=", h)

# for explicit solution
σ = a^2 * τ^2 / h^2
println("σ=", σ)

x = range(0, l, step=h)
print("x:")
println(x)
t = range(0, T, step=τ)
print("t:")
println(t)

mesh = collect(Iterators.product(x, t))
U = sol.(mesh)

er = get_errors(U, U2, N, K)
println("er=", er)

if σ <= 1
    println("Условие Куррента выполнено:", σ, "<= 1\n")
    U2 = explicit(x, t, h, τ)
    srf2 = Plots.surface(x, t, U2', c = :matter)
    plt2 = Plots.plot(x, [U[:, K] U2[:, K]], labels=["точное решение" "явная схема"])
end
srf = PlotlyJS.plot([PlotlyJS.surface(x=x, y=t, z=U, colorscale="Jet"), PlotlyJS.surface(x=x, y=t, z=U2)])

# U3 = implicit_crank(x, t, h, τ, 1)
# plt3 = Plots.plot(x, [U[:, K] U3[:, K]], labels=["точное решение" "неявная схема"])
# U3_crank = implicit_crank(x, t, h, τ, 0.5)
# plt3_crank = Plots.plot(x, [U[:, K] U3_crank[:, K]], labels=["точное решение" "схема Кранка-Николсона"])

# srf = PlotlyJS.plot(PlotlyJS.surface(x=x, y=t, z=U))

# # get errors:
# er1, er2 = get_errors(U, U2, U3, N, K)
# println("ERRORS:")
# println(er1, "\n", er2)

# e1 = [max_abs_error(U[:, i], U2[:, i]) for i in 1:length(t)]
# e2 = [max_abs_error(U[:, i], U3[:, i]) for i in 1:length(t)]
# e3 = [max_abs_error(U[:, i], U3_crank[:, i]) for i in 1:length(t)]
# ep = Plots.plot(t, [e1 e2 e3], labels=["явная" "неявная" "Кранка-Николсона"])
# # ep1 = Plots.plot(t, e1)
# # ep2 = Plots.plot(t, e2)
# # ep3 = Plots.plot(t, e3)
