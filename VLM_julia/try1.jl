
# vortex lattice

# aero
α = 0.1

# planar planform geometry and functions
AR = 1
S = 1
λ = 1
Λ = 0

b = sqrt(AR*S)
c(η) = 2*S/(b*(1+λ)) * (1 .-(1-λ).*abs.(η))
xLE(η) = c(0)/4 .- c(η)./4 .+ tan(Λ)*abs.(η)
xTE(η) = c(0)/4 .+ 3*c(η)./4 .+ tan(Λ)*abs.(η)

xW(ξ, η) = (1 .-ξ).*xLE.(η) + ξ.*xTE(η)
yW(ξ, η) = b/2 * η
zW(ξ, η) = 0*ξ

# M chordwise by N spanwise panels
M = 1;
N = 24;

# spanwise "η" and chordwise "ξ" vortex and control coordinates
ηv = cos.((2*(0:N)) ./ (2*N) * pi)
ηc = cos.((2*(0:N-1) .+ 1) ./ (2*N) * pi)
ξv = sin.((2*(1:M) .- 1) ./ (2*M+1) * pi/2).^2
ξc = sin.((2*(1:M)) ./ (2*M+1) * pi/2).^2

# coordinate vectors to grid matrices
# control points
ηc = reshape(ηc .* ones(N, M), M*N)
ξc = reshape(ξc' .* ones(N, M), M*N)

# left and right bound vortex endpoints
ηvL = reshape(ηv[1:end-1] .* ones(N, M), M*N)
ηvR = reshape(ηv[2:end] .* ones(N, M), M*N)
ξvL = reshape(ξv' .* ones(N, M), M*N)
ξvR = reshape(ξv' .* ones(N, M), M*N)

# planform coordinates to cartesian coordinates
xyzC = [xW(ξc, ηc) yW(ξc, ηc) zW(ξc, ηc)]
xyzVL = [xW(ξvL, ηvL) yW(ξvL, ηvL) zW(ξvL, ηvL)]
xyzVR = [xW(ξvR, ηvR) yW(ξvR, ηvR) zW(ξvR, ηvR)]

n = [zeros(N*M) zeros(N*M) ones(N*M)]
