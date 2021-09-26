using LinearAlgebra

# build vortex lattice geometry from planform parameters
function geomgen(AR, S, λ, Λ, M, N)

    # span, chord function, and LE, TE functions
    b = sqrt(AR*S)
    c(η) = 2*S/(b*(1+λ)) * (1 .-(1-λ).*abs.(η))
    xLE(η) = c(0)/4 .- c(η)./4 .+ tan(Λ)*abs.(b/2*η)
    xTE(η) = c(0)/4 .+ 3*c(η)./4 .+ tan(Λ)*abs.(b/2*η)

    xW(ξ, η) = (1 .-ξ).*xLE.(η) + ξ.*xTE(η)
    yW(ξ, η) = b/2 * η
    zW(ξ, η) = 0*ξ

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

    # left and right trailing edge points
    ξvTE = ones(N*M)

    # planform coordinates to cartesian coordinates
    xyzC = [xW(ξc, ηc) yW(ξc, ηc) zW(ξc, ηc)]
    xyzVL = [xW(ξvL, ηvL) yW(ξvL, ηvL) zW(ξvL, ηvL)]
    xyzVR = [xW(ξvR, ηvR) yW(ξvR, ηvR) zW(ξvR, ηvR)]
    xyzTEL = [xW(ξvTE, ηvL) yW(ξvTE, ηvL) zW(ξvTE, ηvL)]
    xyzTER = [xW(ξvTE, ηvR) yW(ξvTE, ηvR) zW(ξvTE, ηvR)]

    n = [zeros(N*M) zeros(N*M) ones(N*M)]

    return xyzC, xyzVL, xyzVR, xyzTEL, xyzTER, n, c(ηc)

end

################################

# biot savart law on line segment, ray, and line
function biotsavart(r, ra, rb, Rc)
    # unit induced velocity at r by a vortex line from ra to rb
    # bound vortices

    a = r - ra;
    b = r - rb;

    v = 1/(4*pi) * (1/norm(a) + 1/norm(b)) * cross(a, b)./(norm(a)*norm(b) + dot(a, b) .+ Rc.^2)

    return v

end
function sibiotsavart(r, ra, t, Rc)
    # unit induced velocity at r by a semi infinite vortex line from ra along unit vector t
    # trailing vortices
 
    a = r - ra;

    v = 1/(4*pi) * (1/norm(a)) * cross(a, t)./(norm(a) - dot(a, t) .+ Rc.^2)

    return v   
end
function ibiotsavart(r, ra, t, Rc)
    # unit induced velocity at r by an infinite vortex line through ra along unit vector t
    # trefftz plane wake

    a = r - ra;

    v = 1/(2*pi) * cross(a, t)./(dot(a, a) - dot(a, t).^2 .+ Rc.^2)

    return v    
end

################################

# build influence matrix
function vlma(xyzC, xyzVL, xyzVR, xyzTEL, xyzTER, α, β, n)

    Rc = 1e-9
    len = size(n, 1)
    A = zeros(len, len)

    t = [cos(α)*cos(β), -sin(β), sin(α)*cos(β)]
    # slowwwwwwww, but why????
    for i in 1:len
        for j in 1:len
            
            ltv = sibiotsavart(xyzC[i,:], xyzTEL[j,:], t, Rc)
            lbv = biotsavart(xyzC[i,:], xyzTEL[j,:], xyzVL[j,:], Rc)
            bv = biotsavart(xyzC[i,:], xyzVL[j,:], xyzVR[j,:], Rc)
            rbv = biotsavart(xyzC[i,:], xyzVR[j,:], xyzTER[j,:], Rc)
            rtv = sibiotsavart(xyzC[i,:], xyzTER[j,:], t, Rc)

            v = ltv + lbv + bv + rbv - rtv

            A[i, j] = dot(v, n[i,:])
        end
    end
    
    return A
    
end
# build RHS
function vlmb(α, β, n)

    vinf = [cos(α)*cos(β), -sin(β), sin(α)*cos(β)]
    
    len = size(n, 1)
    b = zeros(len)

    for j in 1:len
        b[j] = dot(vinf, n[j,:])
    end

    return b

end

## do calculation as a function

function vlmsolve(α::AbstractVector{T}) where T
    ## remaining variables
    β = T[0.0]
    AR = 6.5
    λ = 0.5
    S = 1.0
    Λ = 0.0
    M = 2
    N = 24

    geom = geomgen(AR, S, λ, Λ, M, N)

    xyzC = T[geom[1]]
    xyzVL = T[geom[2]]
    xyzVR = T[geom[3]]
    xyzTEL = T[geom[4]]
    xyzTER = T[geom[5]]
    n = T[geomT[6]]
    cs = T[geom[7]]

    A = vlma(xyzC, xyzVL, xyzVR, xyzTEL, xyzTER, α, β, n)
    b = vlmb(α, β, n)

    G = A\b
    Gc = vec(sum(reshape(G, N, M), dims=2))


    Cl = 2*Gc./cs[1:N]
    Ccl = Cl.*cs[1:N]*sqrt(AR*S)/S
    CL = sum(2*Gc.*diff(-[xyzVL[1:N,2]; xyzVR[N,2]]))/S

    return CL
end