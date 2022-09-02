function [xiQ, wQ, xiC, wC] = chordLattice(N, LEV)
%chordLattice Generates chordwise discretization, Ltip Rtip true/false
%   Detailed explanation goes here

xw = jacobi(N, 0.5, -0.5, nonzeros(-1*LEV));
xiQ = xw(:,1);
wQ = xw(:,2);

[xiC, wC] = cauchy(xiQ, wQ, 0.5, -0.5);

end

%% Jacobi

function [xw, ab] = jacobi(N, a, b, pn)
%JACOBI Summary of this function goes here
%   Detailed explanation goes here

% N = N-1;

nu=(b-a)/(a+b+2);

if a+b+2 > 128
    mu=exp((a+b+1)*log(2)+((gammaln(a+1)+gammaln(b+1))-gammaln(a+b+2)));
else
	mu=2^(a+b+1)*((gamma(a+1)*gamma(b+1))/gamma(a+b+2));
end



%% recurrence coefficients

n=1:N; nab=2*n+a+b;
A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n=2:N; nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab=[A' [mu; B1; B']];

if N==1
    ab=[nu mu];
end

ab0 = ab;

%% Radau

if length(pn) == 1

    N = N-1;
    if pn == -1 % radau point at -1
      ab(N+1,1)=-1+2*N*(N+a)/((2*N+a+b)*(2*N+a+b+1));
    else % radau point at 1
      ab(N+1,1)=1-2*N*(N+b)/((2*N+a+b)*(2*N+a+b+1));
    end
    N = N+1;
    
end

%% lobatto

if length(pn) == 2

    N = N-2;
    ab(N+2,1)=(a-b)/(2*N+a+b+2);
    ab(N+2,2)=4*(N+a+1)*(N+b+1)*(N+a+b+1)/((2*N+a+b+1)* ...
    (2*N+a+b+2)^2);
    N = N+2;

end
%% Golub-Welsch

J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(ab(n,2));
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];

% w0 = xw(:,2)'

end

%% Cauchy

function [xC, wC] = cauchy(xQ, wQ, a, b)
%CAUCHY Summary of this function goes here
%   Detailed explanation goes here

N = length(xQ);
N2 = N;

%% Cauchy Nodes

% W = @(t) (1-t).^a .* (1+t).^b;

[~, ab] = jacobi(N2, 0, 0, []);

% precomputed principal value of weight function

a = -a; b = -b;
ab2 = a+10*b;

switch ab2
    case -5.5
        C = @(x) 0 + 0.*x;
    case -4.5
        C = @(x) pi + 0.*x;
    case 4.5
        C = @(x) pi + 0.*x;
    case -0.5
        C = @(x) (2*log(sqrt(1-x)./(sqrt(2)+sqrt(1-x)))+log((1+x)./(1-x)))./sqrt(1-x);
    case -5
        C = @(x) (2*log((sqrt(2)+sqrt(1+x))./sqrt(1+x))+log((1+x)./(1-x)))./sqrt(1+x);
    case 0.5
        C = @(x) sqrt(1-x).*(2*log(sqrt(1-x)./(sqrt(2)+sqrt(1-x)))+log((1+x)./(1-x)))+2*sqrt(2);
    case 5
        C = @(x) sqrt(1+x).*(2*log((sqrt(2)+sqrt(1+x))./sqrt(1+x))+log((1+x)./(1-x)))-2*sqrt(2);
    case 0
        C = @(x) log((1+x)./(1-x));
end

wQ = wQ(:)';
xQ = xQ(:);

Q = @(x) wQ*(1./(x-xQ));

f = @(x) C(x) - Q(x);

xC = zeros(N-1, 1);

xQint = [xQ', 1];

for ci = 1:N
    xC(ci) = fzero(f, [xQint(ci)+1e-6, xQint(ci+1)-1e-6]);
end

%% Elhay, Kautsky, Burkardt

wC = cawiq(N2, xC, ones(1,N2), N2, (1:N2), 1, floor((N2+1)/2), ab(1:end-1,1), sqrt(ab(2:end,2)), 0, ab(1,2));

end

%% cawiq
function [ wts, ndx, jdf ] = cawiq ( nt, t, mlt, nwts, ndx, key, nst, aj, bj, jdf, zemu )
  wts = zeros ( nwts, 1 );
  prec = eps;
  if ( nt < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CAWIQ - Fatal error!\n' );
    fprintf ( 1, '  NT < 1.\n' );
    error ( 'CAWIQ - Fatal error!' );
  end
  if ( 1 < nt )
    k = nt - 1;
    for i = 1 : k
      tmp = t(i);
      l = i + 1;
      for j = l : nt
        if ( abs ( tmp - t(j) ) <= prec )
          fprintf ( 1, '\n' );
          fprintf ( 1, 'CAWIQ - Fatal error!\n' );
          fprintf ( 1, '  Knots too close.\n' );
          error ( 'CAWIQ - Fatal error!' );
        end
      end
    end
  end
  l = abs ( key );
  if ( l < 1 || 4 < l )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CAWIQ - Fatal error!\n' );
    fprintf ( 1, '  Magnitude of KEY not between 1 and 4.\n' );
    error ( 'CAWIQ - Fatal error!' );
  end
  k = 1;
  if ( l == 1 )
    for i = 1 : nt
      ndx(i) = k;
      if ( mlt(i) < 1 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'CAWIQ - Fatal error!\n' );
        fprintf ( 1, '  MLT(I) < 1.\n' );
        error ( 'CAWIQ - Fatal error!' );
      end
      k = k + mlt(i);
    end
    n = k - 1;
  elseif ( l == 2 || l == 3 )
    n = 0;
    for i = 1 : nt
      if ( ndx(i) == 0 )
        continue
      end
      if ( mlt(i) < 1 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'CAWIQ - Fatal error!\n' );
        fprintf ( 1, '  MLT(I) < 1.\n' );
        error ( 'CAWIQ - Fatal error!' );
      end
      n = n + mlt(i);
      if ( ndx(i) < 0 && l == 3 )
        continue
      end
      ndx(i) = abs ( k ) * i4_sign ( ndx(i) );
      k = k + mlt(i);
    end
    if ( nwts + 1 < k )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'CAWIQ - Fatal error!\n' );
      fprintf ( 1, '  NWTS + 1 < K.\n' );
      error ( 'CAWIQ - Fatal error!' );
    end
  elseif ( l == 4 )
    for i = 1 : nt
      ip = abs ( ndx(i) );
      if ( ip == 0 )
        continue
      end
      ipm = ip + mlt(i);
      if ( nwts < ipm )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'CAWIQ - Fatal error!\n' );
        fprintf ( 1, '  NWTS < IPM.\n' );
        error ( 'CAWIQ - Fatal error!' );
      end
      if ( i == nt )
        break
      end
      l = i + 1;
      for j = l : nt
        jp = abs ( ndx(j) );
        if ( jp ~= 0 )
          if ( jp <= ipm && ip <= jp + mlt(j) )
            break
          end
        end
      end
    end
  end
  if ( nst < floor ( ( n + 1 ) / 2 ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CAWIQ - Fatal error!\n' );
    fprintf ( 1, '  NST < ( N + 1 ) / 2.\n' );
    error ( 'CAWIQ - Fatal error!' );
  end
  if ( zemu <= 0.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CAWIQ - Fatal error!\n' );
    fprintf ( 1, '  ZEMU <= 0.\n' );
    error ( 'CAWIQ - Fatal error!' );
  end
  if ( n <= 1 )
    for i = 1 : nt
      if ( 0 < ndx(i) )
        wts( abs ( ndx(i) ) ) = zemu;
        return
      end
    end
  end
  if ( jdf == 0 )
    z = zeros(1,nst);
    z(1) = 1.0;
    [ aj, z ] = imtqlx ( nst, aj, bj, z );
    jdf =1;
    for i = 1 : nst
      bj(i) = z(i) * z(i);
    end
  end
  for i = 1 : nt
    if ( ndx(i) <= 0 )
      cycle
    end
    m = mlt(i);
    nm = n - m;
    mnm = max ( nm, 1 );
    l = min ( m, nm + 1 );
    xk = zeros(mnm,1);
    k = 1;
    for j = 1 : nt
      if ( ndx(j) ~= 0 )
        if ( j ~= i )
          for jj = 1 : mlt(j)
            xk(k) = t(j);
            k = k + 1;
          end
        end
      end
    end
    r = zeros(l,1);
    r(1) = 1.0 / zemu;
    k = ndx(i);
    wts(k:k+m-1) = cwiqd ( m, mnm, l, t(i), xk, nst, aj, bj, r );
    if ( key < 0 )
      continue
    end
    tmp = 1.0;
    for j = 2 : m - 1
      p = j;
      tmp = tmp * p;
      l = k + j;
      wts(l) = wts(l) / tmp;
    end
  end
  return
end

%% cwiqd
function d = cwiqd ( m, nm, l, v, xk, nstar, phi, a, r )
  wf = zeros ( nstar, 1 );
  for j = 1 : nstar
    wf(j) = a(j);
    for i = 1 : nm
      wf(j) = wf(j) * ( phi(j) - xk(i) );
    end
  end
  y = zeros ( m );
  for i = 1 : m
    sum = 0.0;
    for j = 1 : nstar
      sum = sum + wf(j);
      wf(j) = wf(j) * ( phi(j) - v );
    end
    y(i) = sum;
  end
  for i = 1 : nm

    tmp = v - xk(i);

    last = min ( l, i + 1 );
    for jr = 2 : last
      j = last - jr + 2;
      r(j) = tmp * r(j) + r(j-1);
    end

    r(1) = tmp * r(1);

  end
  d(m) = y(m) / r(1);

  if ( m == 1 )
    return
  end
  z = zeros ( m );
  z(1) = 1.0 / r(1);
  for i = 2 : m
    sum = 0.0;
    minil = min ( i, l );
    for j = 2 : minil
      k = i - j + 1;
      sum = sum + r(j) * z(k);
    end
    z(i) = - sum * z(1);
  end
  for i = 2 : m
    sum = 0.0;
    for j = 1 : i
      k = m - i + j;
      sum = sum + z(j) * y(k);
    end
    k = m - i + 1;
    d(k) = sum;
  end
  return
end

%% imtqlx
function [ d, z ] = imtqlx ( n, d, e, z )
  itn = 30;
  prec = eps;
  if ( n == 1 )
    return
  end
  e(n) = 0.0;
  for l = 1 : n
    j = 0;
    while ( 1 )
      for m = l : n
        if ( m == n )
          break
        end
        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) )
          break
        end
      end
      p = d(l);
      if ( m == l )
        break
      end
      if ( j == itn )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'IMTQLX - Fatal error!\n' );
        fprintf ( 1, '  Iteration limit exceeded.\n' );
        error ( 'IMTQLX - Fatal error!' );
      end
      j = j + 1;
      g = ( d(l+1) - p ) / ( 2.0 * e(l) );
      r =  sqrt ( g * g + 1.0 );
      g = d(m) - p + e(l) / ( g + r8_sign ( g ) * abs ( r ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;
      for ii = 1 : mml
        i = m - ii;
        f = s * e(i);
        b = c * e(i);
        if ( abs ( f ) >= abs ( g ) )
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e(i+1) = f * r;
          s = 1.0 / r;
          c = c * s;
        else
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e(i+1) = g * r;
          c = 1.0 / r;
          s = s * c;
        end
        g = d(i+1) - p;
        r = ( d(i) - g ) * s + 2.0 * c * b;
        p = s * r;
        d(i+1) = g + p;
        g = c * r - b;
        f = z(i+1);
        z(i+1) = s * z(i) + c * f;
        z(i) = c * z(i) - s * f;
      end
      d(l) = d(l) - p;
      e(l) = g;
      e(m) = 0.0;
    end
  end
  for ii = 2 : n
     i = ii - 1;
     k = i;
     p = d(i);
     for j = ii : n
       if ( d(j) < p )
         k = j;
         p = d(j);
       end
     end
     if ( k ~= i )
       d(k) = d(i);
       d(i) = p;
       p = z(i);
       z(i) = z(k);
       z(k) = p;
     end
  end
  return
end

%% r8_sign
function value = r8_sign ( x )
  if ( 0 <= x )
    value = +1.0;
  else
    value = -1.0;
  end
  return
end
