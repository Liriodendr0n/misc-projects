clear
close all
clc

%% velocities

syms x y z Gamma U

% semi-infinite doublet line
phiDL = Gamma/(4*pi) * z/(y^2+z^2) * x/sqrt(x^2+y^2+z^2);
uvwDL = gradient(phiDL, [x; y; z]);

% vorton (to investigate lift contributions???)
% all the lift *should* come from the "bound end" of the doublet line
AV = [0; Gamma/(4*pi) * 1/sqrt(x^2+y^2+z^2); 0];
uvwV = curl(AV, [x; y; z]);

uvw = simplify(uvwDL - 0*uvwV);

%% pressure disturbances

p = -2*U*uvw(1) - uvw(1)^2 - uvw(2)^2 - uvw(3)^2;

%% control volume boundaries

X = 1;
Y = 1;
Z = 1;

% relevant flow scales
G = 0.1;
Uinf = 10;

%% vertical momentum integration
% front face x = -X
m_front = -Uinf*vpaintegral(vpaintegral(subs(uvw(3), [x, Gamma], [-X, G]), y, [-Y, Y], 'abstol', 1e-6), z, [-Z, Z], 'abstol', 1e-6);
% back face x = X
m_back = -Uinf*vpaintegral(vpaintegral(subs(uvw(3), [x, Gamma], [X, G]), y, [-Y, Y], 'abstol', 1e-6), z, [-Z, Z], 'abstol', 1e-6);

%% vertical pressure integration
% top face z = Z
p_top = 0.5*vpaintegral(vpaintegral(subs(p, [Gamma, U, z], [G, Uinf, Z]), y, [-Y, Y], 'abstol', 1e-6), x, [-X, X], 'abstol', 1e-6);
% bottom face z = -Z
p_bottom = 0.5*vpaintegral(vpaintegral(subs(p, [Gamma, U, z], [G, Uinf, -Z]), y, [-Y, Y], 'abstol', 1e-6), x, [-X, X], 'abstol', 1e-6);

%% lift contributions
Lm = round(double((m_back - m_front)), 6)
Lp = round(double((p_bottom - p_top)), 6)
L = Lm + Lp
LKJ = Uinf*G;

%%
