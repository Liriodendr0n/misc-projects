clear
close all
clc

alpha = deg2rad(5);
xyz = [0; 0; 0];
N = 301;

S = 1;
AR = 1;
Gamma = deg2rad(0);

lambdas = linspace(0, 1, 30);
Lambdas = deg2rad(0);

for l = 1:length(lambdas)
    parfor L = 1:length(Lambdas)
        
        lambda = lambdas(l);
        Lambda = Lambdas(L);
        
        b = sqrt(AR*S);
        c = S/b;
        %% build geometry
        [xyzV, xyzC, xyzD, cCs, cVs, ns] = wingGeom(xyz, N, AR, S, lambda, Lambda, Gamma);

        %% build and solve linear system
        RHS = weisslRHS(alpha, ns, 0*ones(N,1));

        A = weisslA(xyzC, xyzV, ns);
        Atrefftz = weisslAtrefftz(xyzD, xyzV, ns);

        G = A\RHS;

        % local Cl
        Cl = 2*G./cCs;
        % local Loading
        Ccl = 2*G/c;
        % perpendicular Cl
        Clp = Cl/(cos(Lambda)^2);

        %% downwash
        % local induced angle (includes bound vortex)
        alphai = alpha - Cl/(2*pi);

        % trefftz downwash
        w = -0.5*Atrefftz*G;

        Cdi = w.*Cl;
        Ccdi = w.*Ccl;

        %% results
        CL = trapz(xyzC(2,:)/b, Ccl);
        CDi = trapz(xyzC(2,:)/b, Ccdi);

        e(L,l) = (CL^2)/(pi*AR*CDi);
    end
end


figure
hold on
for i = 1:length(Lambdas)
    plot(lambdas, e(i,:))
end
[maxv, maxi] = max(e(:,:));
lambdas(maxi)
        