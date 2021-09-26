using Plots
using UnicodePlots

#include("vlmUtils.jl")


α = deg2rad(10)
β = deg2rad(0)

AR = 6.5
S = 1
λ = 0.5
Λ = deg2rad(0)

# M chordwise by N spanwise (right now only 1 chordwise panel)
M = 2
N = 24

geom = geomgen(AR, S, λ, Λ, M, N)

xyzC = geom[1]
xyzVL = geom[2]
xyzVR = geom[3]
xyzTEL = geom[4]
xyzTER = geom[5]
n = geom[6]
cs = geom[7]

@time begin
A = vlma(xyzC, xyzVL, xyzVR, xyzTEL, xyzTER, α, β, n)
b = vlmb(α, β, n)

G = A\b
Gc = vec(sum(reshape(G, N, M), dims=2))


Cl = 2*Gc./cs[1:N]
Ccl = Cl.*cs[1:N]*sqrt(AR*S)/S
CL = sum(2*Gc.*diff(-[xyzVL[1:N,2]; xyzVR[N,2]]))/S



end
##
##

#plot(2*xyzC[1:N,2]/(sqrt(AR*S)), Cl)
plt = lineplot(2*xyzC[1:N,2]/(sqrt(AR*S)), Cl, ylim = [0 round(maximum([Cl; Ccl])*1.25, digits = 1)], xlabel = "η", ylabel = "Cl, Ccl", name = "Cl")
lineplot!(plt, 2*xyzC[1:N,2]/(sqrt(AR*S)), Ccl, name = "Ccl")