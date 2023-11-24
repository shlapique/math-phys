using Plots
using Unzip

# u(x, 0) = sin(x)
function psi(x)
    return sin(x)
end

# u(0, t)
function phi0(t)
    return exp(-0.5 * t)
end

# u(pi, t)
function phi1(t)
    return -exp(-0.5 * t)
end

# U(x, t) = exp(-0,5t)sin(x)
function sol(x, t)
    return exp(-0.5 * t) * sin.(x)
end

function run_through(A, b)
    n = size(A)[1]
    print(n)
    x = zeros(n)
    P = zeros(n)
    Q = zeros(n)
    P[1] = -A[1, 2] / A[1, 1]
    Q[1] = b[1] / A[1, 1]
    for i = 2:n
        if i == n
            P[i] = 0
        else
            println("i=", i)
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

function explicit(K, t, τ, σ, h, x, eps)
    # N = length(x)
    U = zeros(K, N)

    for j in 1:N
        # inital cond
        U[1, j] = psi(x[j])
    end
    for k in 2:length(t)-1
        U[k, 1] = phi0(t[k])
        for j in 2:length(x)-1
            U[k, j] = σ * U[k-1, j-1] + (1 - 2*σ) * U[k-1, j] + σ * U[k-1, j+1]
        end
        U[k, end] = phi1(t[k])
        # if eps == 1
        #     U[k + 1, 0] = U[k + 1, 1] - h * phi0(t)
        #     U[k + 1, N - 1] = phi0(t)
        # elseif eps == 2
        #     U[k + 1, 0] = (2 * h * phi0(t) - 4 * U[k + 1, 1] + U[k + 1, 2]) / (2 * h - 3)
        #     U[k + 1, N - 1] = phi0(t)
        # elseif eps == 3
        #     U[k + 1, 0] = (phi0(t) * h * τ * 2 - U[k + 1, 1] * (2 * τ) - U[k, 0] * h^2) / (-2 * τ - h^2 + c * τ * h^2 + h * τ * 2)
        #     U[k + 1, N - 1] = phi0(t)
        # end
    end
    return U
end
            

l = pi
a = 1

N = 50
K = 7000
T = 3 # sec
τ = T / K

h = (l - 0) / N

# for explicit solution
σ = a^2 * τ / h^2

println("σ=", σ)

# get analytical solution
x = range(0, l, step=h)
t = range(0, T, step=τ)

# println("Enter t:")
# t_cur = readline()
# t_cur = parse(Float64, t_cur)
t_cur = 0 
U = sol(x, t_cur)
display(U)

if a * τ / h^2 <= 0.5
    println("Условие Куррента выполнено:", a * τ / h^2, "<= 0.5\n")
    U2 = explicit(K, t, τ, σ, h, x, eps)
    plot(x, U2[3000, :])
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
