using Plots

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
    return exp(-0.5 * t)
end

# U(x, t) = exp(-0,5t)sin(x)
function sol(x, t)
    return exp(-0.5 * t) * sin(x)
end

function solve_sol(l, r, eps)
    x_range = range(0, stop=2pi, length=10000)
end

# function run_through(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64},
#                 d::Vector{Float64}, n::Int64)
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

x_range = range(0, stop=2pi, length=10000)

# test
A = [-10 9 0 0 0
    -5 -21 -8 0 0
    0 7 12 2 0
    0 0 0 8 2
    0 0 0 2 10]
display(A)
b = [7 29 31 56 -24]
println(b)
test = run_through(A, b)
println("test\n", test)
