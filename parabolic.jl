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
function sol(x, t)
    return exp(-0.5 * t) * sin.(x)
end

function f(x, t)
    return 0.5 * exp(-0.5 * t) * sin(x)
end


function run_through(A, b)
    n = size(A)[1]
    # print(n)
    x = zeros(n)
    P = zeros(n)
    Q = zeros(n)
    P[1] = -A[1, 2] / A[1, 1]
    Q[1] = b[1] / A[1, 1]
    for i = 2:n
        if i == n
            P[i] = 0
        else
            # println("i=", i)
            P[i] = -A[i, i + 1] / (A[i, i] + A[i, i - 1] * P[i - 1])
        end
        Q[i] = (b[i] - A[i, i - 1] * Q[i - 1]) / (A[i, i] + A[i, i - 1] * P[i - 1])
    end
    x[n] = Q[n]

    for i = n-1:-1:1
        x[i] = P[i] * x[i + 1] + Q[i]
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

function implicit(x, t, h, τ, θ)
    U = zeros(length(x), length(t))
    σ = τ / h^2
    # inital cond
    for i in 1:length(x)
        U[i, 1] = ψ(x[i])
    end
	for i in 2:length(t)
			A = zeros(length(x)-2, length(x)-2)
            n = size(A)[1]
			A[1, 1] = -(1 + 2*σ)
			A[1, 2] = σ
            # for j in 2:length(A)-1
            for j in 2:n-1
                A[j, j-1] = σ
                A[j, j] = -(1 + 2*σ)
                A[j, j+1] = σ
            end
            A[end, end-1] = σ
            A[end, end] = -(1 + 2*σ)

            b = -U[2:end-1, i-1]

			b[1] -= σ * ϕ0(t[i])
            b[end] -= σ * ϕ1(t[i])

			U[1, i] = ϕ0(t[i])
            U[end, i] = ϕ1(t[i])
            U[2:end-1, i] = run_through(A, b)
	end
	return U
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
U = sol(x, t_cur)
display(U)

if a * τ / h^2 <= 0.5
    println("Условие Куррента выполнено:", a * τ / h^2, "<= 0.5\n")
    U2 = explicit(x, t, h, τ)
    U3 = implicit(x, t, h, τ, 1)
    plot(x, U2[:, K])
    plot!(x, U3[:, K])
end

plot!(x, U)



# mesh = collect(Iterators.product(x, t))
# display(mesh)

# t1, t2 = unzip(mesh)
# t1 = t1[:,1]
# t2 = t2[1,:]
# display(t1)
# display(t2)
# plot(t1, t2)

# y = sol.(mesh)
# display(uz)
# typeof(uz)
# display(y)
# plot(x, y)
