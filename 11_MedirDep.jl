### Medidas de dependencia y concordancia
### ⟶ Continuación de 09_Simular.jl
### Autor: Arturo Erdely
### Fecha: 2026-01-29


using Plots, LaTeXStrings, Distributions
using HCubature # para integración numérica de funciones de varias variables


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
    n = 10_000
    U, V = simular_copula(n, θ)
    scatter(U, V, xlabel = L"U", ylabel = L"V", legend = false, size = (500, 500), markersize = 0.5)
end

# Cópula empírica
Cn(u,v) = (1/n) * sum( (U .≤ u) .* (V .≤ v) )


## Medidas de dependencia (Schweizer-Wolff) y concordancia (Spearman)

# Teóricas 
function σ(θ) # Schweizer-Wolff
    integrando(uv) = abs(C(uv[1], uv[2], θ) - uv[1]*uv[2])
    valor = hcubature(integrando, [0.0, 0.0], [1.0, 1.0])[1]
    return 12 * valor
end
function ρ(θ) # Spearman
    integrando(uv) = C(uv[1], uv[2], θ) - uv[1]*uv[2]
    valor = hcubature(integrando, [0.0, 0.0], [1.0, 1.0])[1]
    return 12 * valor
end

# Empíricas
function σn(m = n) # Schweizer-Wolff empírico
    suma = 0.0
    u = collect(range(0.0, 1.0, length = m))
    v = collect(range(0.0, 1.0, length = m))
    for i ∈ 1:m
        for j ∈ 1:m
            suma += abs(Cn(u[i], v[j]) - u[i]*v[j])
        end
    end
    return (12/(m^2)) * suma
end
function ρn(m = n) # Spearman empírico
    suma = 0.0
    u = collect(range(0.0, 1.0, length = m))
    v = collect(range(0.0, 1.0, length = m))
    for i ∈ 1:m
        for j ∈ 1:m
            suma += Cn(u[i], v[j]) - u[i]*v[j]
        end
    end
    return (12/(m^2)) * suma
end


## Cálculo de las medidas teóricas y empíricas
@show θ;
@show σ(θ);
@show ρ(θ);
@show σn(1_000);
@show ρn(1_000);
