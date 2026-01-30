### Estimación paramétrica de una cópula
### cuando no existe densidad conjunta
### ⟶ Continuación de 09_Simular.jl
### Autor: Arturo Erdely
### Fecha: 2026-01-29


using Plots, LaTeXStrings, Distributions
using Optim # para optimización numérica


## Cópula arquimediana de 07_CopArquim.jl

# generador
ϕ(t,θ) = (1.0 - t) ^ θ # 0 ≤ t ≤ 1, θ ≥ 1

# pseudo-inversa del generador
function ϕi(z, θ) # 0 ≤ z ≤ ∞
    if 0.0 ≤ z ≤ 1.0 
        return 1.0 - z^(1/θ)
    else
        return 0.0
    end
end

# cópula arquimediana
C(u,v,θ) = ϕi(ϕ(u,θ) + ϕ(v,θ),θ)

# curva frontera de la región cero 
g(u,θ) = 1 - (1 - (1-u)^θ)^(1/θ) 


## Simular cópula arquimediana de 08_Simular.jl

# ψ⁻¹(t|u) como función de t, para u y θ fijos
function ψi(t, u, θ)
    if t > (1-u)^(θ-1) 
        return 1 - (1-u)*(1/(t^(θ/(θ-1))) - 1)^(1/θ)
    else
        return g(u, θ)
    end
end

# Función para simular n observaciones de la cópula arquimediana con parámetro θ
function simular_copula(n, θ)
    U = rand(n)
    V = similar(U)
    for i ∈ 1:n
        u = U[i]
        w = rand()
        V[i] = ψi(w, u, θ)
    end
    return U, V
end

# Simulación de n observaciones de la cópula arquimediana con parámetro θ
begin
    θ = 2.3
    n = 1_000
    U, V = simular_copula(n, θ)
    scatter(U, V, xlabel = L"U", ylabel = L"V", legend = false, size = (500, 500), markersize = 0.5)
end

# Cópula empírica
Cn(u,v) = (1/n) * sum( (U .≤ u) .* (V .≤ v) )


## Enfoque clásico: Mínimos cuadrados ordinarios 

# Función objetivo para MCO
function obj_mco(θ, uu, vv)
    suma = 0.0
    for i ∈ eachindex(uu)
        suma += (Cn(uu[i], vv[i]) - C(uu[i], vv[i], θ))^2
    end
    return suma
end

obj_mco(1.0, U, V)
obj_mco(2.3, U, V)
obj_mco(5.0, U, V)

# Optimización numérica para obtener el estimador MCO de θ

@doc Optim

# https://julianlsolvers.github.io/Optim.jl/stable/ 

h(θ) = obj_mco(θ[1], U, V)
res_mco = optimize(h, [1.0])
# estimación puntual
θ_mco = Optim.minimizer(res_mco)


## Enfoque bayesiano: Approximate Bayesian Computation (ABC)

@time begin # 15 seg aprox con m = 10_000
    m = 10_000
    k = 100
    n = length(U)
    Cn_vec = zeros(n)
    for i ∈ 1:n 
        Cn_vec[i] = Cn(U[i], V[i]) # precálculo de la cópula empírica
    end
    θsim = rand(Uniform(1.0, 10.0), m) # Afinar el intervalo según el caso
    ε = zeros(m)
    for j ∈ 1:m
        Usim, Vsim = simular_copula(n, θsim[j])
        CnSim(u,v) = (1/n) * sum( (Usim .≤ u) .* (Vsim .≤ v) )
        dist = 0.0
        for i ∈ 1:n
            dist += (Cn_vec[i] - CnSim(U[i], V[i]))^2
        end
        ε[j] = dist
    end
    idx = sortperm(ε)[1:k]
    θ_abc = θsim[idx]
end;

# estimación puntual
mean(θ_abc)

# estimación de intervalo de probabilidad al 95%
quantile(θ_abc, [0.025, 0.975])

# valores extremos
extrema(θ_abc)

# histograma de la densidad a posteriori
histogram(θ_abc, bins = Int(round(√k)), xlabel = L"\theta", ylabel = "Frecuencia", legend = false, size = (500, 400))


## Rotaciones 

# original (U, V) ~ C(u,v)
scatter(U, V, xlabel = L"U", ylabel = L"V", legend = false, size = (500, 500), markersize = 0.5)

# (1-U, V) ~ C*(u,v) = u - C(u, 1-v)
scatter(1 .- U, V, xlabel = L"1 - U", ylabel = L"V", legend = false, size = (500, 500), markersize = 0.5)

# (1-U, 1-V) ~ C**(u,v) = u + v - 1 + C(1-u, 1-v) 
scatter(1 .- U, 1 .- V, xlabel = L"1 - U", ylabel = L"1 - V", legend = false, size = (500, 500), markersize = 0.5)

# (U, 1-V) ~ C***(u,v) = v - C(1-u, v)
scatter(U, 1 .- V, xlabel = L"U", ylabel = L"1 - V", legend = false, size = (500, 500), markersize = 0.5)
