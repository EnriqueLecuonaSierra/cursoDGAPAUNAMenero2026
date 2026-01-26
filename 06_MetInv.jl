### Método de inversión: Extraer cópula subyacente
### Autor: Arturo Erdely
### Fecha: 2026-01-25


## Paquetes instalados previamente

using Plots, LaTeXStrings, Distributions
using HCubature # integración numérica en varias dimensiones



## Probando el paquete HCubature

@doc HCubature

# ∫∫xydxdy = 1/4 sobre [0,1]x[0,1]
g(x) = x[1]*x[2]
hcubature(g, [0.0, 0.0], [1.0, 1.0]) # integración de g sobre el cuadrado unitario



## Función auxiliar: evaluar función sobre malla 
function mallaXY(fXY, xmin, xmax, ymin, ymax, nx, ny = nx)
    x = collect(range(xmin, xmax, length = nx))
    y = collect(range(ymin, ymax, length = ny))
    matriz = zeros(ny, nx)
    for i ∈ 1:nx, j ∈ 1:ny
        matriz[j, i] = fXY(x[i], y[j])
    end
    return (X = x, Y = y, Z = matriz)
end



## Extrayendo la cópula Gaussiana

@doc MvNormal 

function dcópula_gaussiana(r) # densidad de la cópula gaussiana
    Σ = [1.0 r; r 1.0]
    N2 = MvNormal([0.0, 0.0], Σ)
    Z = Normal(0, 1)
    c(u,v) = pdf(N2, [quantile(Z, u), quantile(Z, v)]) / (pdf(Z, quantile(Z, u)) * pdf(Z, quantile(Z, v)))
    return c
end

function cópula_gaussiana(r)
    c = dcópula_gaussiana(r)
    cc(x) = c(x[1], x[2])
    C(u,v) = hcubature(cc, [0.0, 0.0], [u, v])[1]
    return C
end


r = -0.8
C = cópula_gaussiana(r)
c = dcópula_gaussiana(r)
@time mC = mallaXY(C, 0, 1, 0, 1, 30);
mc = mallaXY(c, 0.01, 0.99, 0.01, 0.99, 100);
sC = surface(mC.X, mC.Y, mC.Z, xlabel = L"u", ylabel = L"v", zlabel = L"C(u,v)", 
             title = "Cópula Gaussiana (r = $r)")
sc = surface(mc.X, mc.Y, mc.Z, xlabel = L"u", ylabel = L"v", zlabel = L"c(u,v)", 
             title = "Densidad cópula Gaussiana (r = $r)")
nivC = contour(mC.X, mC.Y, mC.Z, xlabel = L"u", ylabel = L"v", legend = false,
               title = "Cópula Gaussiana (r = $r)", fill = true, size=(420,400))
nivc = contour(mc.X, mc.Y, mc.Z, xlabel = L"u", ylabel = L"v", 
               title = "Densidad cópula Gaussiana (r = $r)", fill = true, size=(420,400))
#

# Comparar versus W, Π y M 

begin
    n = 100
    W(u,v) = max(u + v - 1, 0) # cópula cota inferior FH
    Π(u,v) = u * v             # cópula independencia
    M(u,v) = min(u, v)         # cópula cota superior FH
    mW = mallaXY(W, 0, 1, 0, 1, n)
    mΠ = mallaXY(Π, 0, 1, 0, 1, n)
    mM = mallaXY(M, 0, 1, 0, 1, n)
    nivW = contour(mW.X, mW.Y, mW.Z, xlabel = L"u", ylabel = L"v", size = (420, 400), 
                   fill = true, title = L"W(u,v)", legend = false)
    nivΠ = contour(mΠ.X, mΠ.Y, mΠ.Z, xlabel = L"u", ylabel = L"v", size = (420, 400), 
                   fill = true, title = L"\Pi(u,v)", legend = false)
    nivM = contour(mM.X, mM.Y, mM.Z, xlabel = L"u", ylabel = L"v", size = (420, 400), 
                   fill = true, title = L"M(u,v)", legend = false)
end;

@time begin
    r = 0.8 # probar r ∈ {-0.8, 0.8, 0.0}
    C = cópula_gaussiana(r)
    mC = mallaXY(C, 0, 1, 0, 1, 30);
    nivC = contour(mC.X, mC.Y, mC.Z, xlabel = L"u", ylabel = L"v", legend = false,
                   title = "Cópula Gaussiana (r = $r)", fill = true, size=(420,400))
end;

nivC
nivW
nivΠ
nivM



## Cópula Gaussiana + marginales arbitrarias 

#=

    (X,Y) ~ f(x,y) = c(F_X(x), F_Y(y))⋅f_X(x)⋅f_Y(y)

    f = densidad conjunta
    c = densidad de la cópula
    F_X, F_Y = funciones de distribución marginales
    f_X, f_Y = densidades marginales
=#

function dSklar(c, X, Y)
    f(x,y) = c(cdf(X, x), cdf(Y, y)) * pdf(X, x) * pdf(Y, y)
    return f
end


# Marginales: Beta(5,2) y Exponencial(2)
begin
    r = 0.7
    c = dcópula_gaussiana(r)
    X = Beta(5, 2)
    Y = Exponential(2)
    f = dSklar(c, X, Y)
    mf = mallaXY(f, 0.0, 1.0, 0.0, 6.0, 50)
end;
sf = surface(mf.X, mf.Y, mf.Z, xlabel = L"x", ylabel = L"y", zlabel = L"f(x,y)", 
                 title = "Densidad conjunta (r = $r)")
cf = contour(mf.X, mf.Y, mf.Z, xlabel = L"x", ylabel = L"y", 
                 title = "Densidad conjunta (r = $r)", fill = true, size=(420,400))
#
ff(x) = f(x[1],x[2])
hcubature(ff, [0, 0], [1, 10]) # debe ser cercano a 1 


# Marginales: Normales(0,1) ⟶ ¡Normal bivariada!
begin
    r = 0.7
    c = dcópula_gaussiana(r)
    X = Normal(0, 1)
    Y = Normal(0, 1)
    f = dSklar(c, X, Y)
    mf = mallaXY(f, -3, 3, -3, 3, 50)
end;
sf = surface(mf.X, mf.Y, mf.Z, xlabel = L"x", ylabel = L"y", zlabel = L"f(x,y)", 
                 title = "Densidad conjunta (r = $r)")
cf = contour(mf.X, mf.Y, mf.Z, xlabel = L"x", ylabel = L"y", 
                 title = "Densidad conjunta (r = $r)", fill = true, size=(420,400))
#
