clear
close all
clc

syms z zeta zeta0 zeta1

%%

h = 1;
tau = 1i;

%% free flow

Gsf(z, zeta) = 1/(sym(2*pi)*(z-zeta));
Ggf = 1i*Gsf;

%% cascade

Hsc(z, zeta) = 1/(2*tau) * (cot(pi*(z-zeta)/tau) - tau/(pi*(z-zeta))) + 1i/(2*tau);
Hgc = 1i*Hsc;

Gsc = Gsf + Hsc;
Ggc = Ggf + Hgc;

%% tunnel

B(z, zeta) = 1/(2*h) * (exp(pi*(z-zeta)/h) - 1)^(-1) - 1/(2*pi*(z-zeta));
E(z, zeta) = -1/(2*h) * (exp(pi*(z-conj(zeta))/h) + 1)^(-1);

Hst = B + E + 1/h;
Hgt = 1i*(B - E);

Gst = Gsf + Hst;
Ggt = Ggf + Hgt;

%%

syms x y

figure
hold on
fsurf(angle(Ggc(x+1i*y, 0)), 'edgecolor', 'none')

