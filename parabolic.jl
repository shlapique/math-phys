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

function get_sol(x, t)
    solution = zeros(0)
    for i in x
        for j in t
            append!(solution, sol(i, j))
        end
    end
    return solution
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
    println("KROL")

    for i = n-1:-1:1
        x[i] = P[i] * x[i + 1] + Q[i]
    end
    return x
end

function explicit(K, t, τ, h, a, c, x, eps)

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
t_cur = 1 
U = sol(x, t_cur)
display(U)

plot(x, U)

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
