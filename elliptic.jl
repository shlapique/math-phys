using PlotlyJS
using Plots
using Unzip

function ϕ1(y)
    return cos(y)
end


function ϕ2(y)
    return 0
end


function ϕ3(x)
    return cos(x)
end


function ϕ4(x)
    return 0
end


function sol(data)
    return cos(data[1]) * cos(data[2])
end

function liebmann(x, y, hx, hy, ε)
    K = length(x)
    N = length(y)

    U = zeros(length(x), length(y))
    # lets find ∂Ω
    for i in 1:N
        U[1, i] = ϕ1(y[i])
        U[end, i] = ϕ2(y[i])
    end
    for j in 1:K
        U[j, 1] = ϕ3(x[j])
        U[j, end] = ϕ4(x[j])
    end

    count = 0
    max_diff = ε + 1.0
    while max_diff > ε
        max_diff = 0
        for i in 2:K-1
            for j in 2:N-1
                new_val = (
                           (U[i + 1, j] * (hy^2) +
                             U[i - 1, j] * (hy^2) +
                             U[i, j + 1] * (hx^2) +
                             U[i, j - 1] * (hx^2)) / (2*hy^2 + 2*hx^2 - 2*hx^2 * hy^2)
                          )
                diff = abs(new_val - U[i, j])
                max_diff = max(max_diff, diff)
                U[i, j] = new_val
            end
        end
        count += 1
    end
    return U, count
end

function get_errors(U, U2, K, N)
    err = 0
    for i in 1:K
        for j in 1:N
            err += (U[i, j] - U2[i, j])^2
        end
    end
    return err/((K+1)*(N+1))
end

function max_abs_error(A, B)
    return maximum(abs.(A - B))
end


lx = pi/2
ly = pi/2
K = 20 
N = 20

X = range(0, ly, length=N)

hx = lx / (K-1)
hy = ly / (N-1)

x = range(0, lx, step=hx)
y = range(0, ly, step=hy)

ε = 1e-5


mesh = collect(Iterators.product(x, y))
U = sol.(mesh)
# srf = Plots.surface(x, y, U, c = :matter)
srf = PlotlyJS.plot(PlotlyJS.surface(x=x, y=y, z=U, colorscale="Jet"))

U2, count = liebmann(x, y, hx, hy, ε)
srf2 = PlotlyJS.plot(PlotlyJS.surface(x=x, y=y, z=U2, colorscale="Blackbody"))

srfu = PlotlyJS.plot([PlotlyJS.surface(x=x, y=y, z=U2, colorscale="Blackbody"),
                      PlotlyJS.surface(x=x, y=y, z=U, colorscale="Jet")])

er1 = get_errors(U, U2, K, N)

e = [max_abs_error(U[:, i], U2[:, i]) for i in 1:length(y)]
ep = Plots.plot(x, e, title="график погрешности от h_y")
e_x = [max_abs_error(U[:, i], U2[:, i]) for i in 1:length(x)]
ep_x = Plots.plot(x, e_x, title="график погрешности от h_x")
