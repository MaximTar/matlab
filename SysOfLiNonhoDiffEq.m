close all
clear
clc

%% initialization
syms th(t) f(t) r(t) u v % th - theta, f - phi
assume(t, 'real')

R = 6700000;
mu = 398600.4415*1E9;

ma = 1; % ma << mb
mb = 10000;
mpr = mb * ma / (mb + ma);
om = sqrt(mu / R^3);

Fx = 0;
Fy = 0;
Fz = 0;

F = [cos(th)*cos(f) sin(th)*cos(f) sin(f); -sin(th) cos(th) 0; -cos(th)*sin(f) -sin(th)*sin(f) cos(f)] * [Fx; Fy; Fz];

Fmat = formula(F);

Fr  = Fmat(1,:);
Fth = Fmat(2,:);
Ff  = Fmat(3,:);

%syms ma om 
Qth = -Fth/(ma*om^2*r*cos(f));
Qf  = -Ff /(ma*om^2*r);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
Qr  = -Fr /(ma*om^2);

syms T ma om r0 dr0

lmbd = T/(ma*om^2*r);

ode1 = diff(th,2) + 2 * (diff(th) + 1) * (diff(r) / r - diff(f) * tan(f)) + 3 * sin(th) * cos(th) == Qth;
ode2 = diff(f,2) + 2 * diff(f) * (diff(r) / r) + sin(f) * cos(f) * ((diff(th) + 1)^2 + 3 * cos(th)^2) == Qf;
ode3 = diff(r,2) + r * (lmbd - (diff(f)^2 + (diff(th) + 1)^2 * cos(f)^2 + 3 * cos(f)^2 * cos(th)^2 - 1))== Qr;
odes = [ode1; ode2; ode3];

%% order reduction
syms x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) xth2 xf2 xr2
x = [x1(t) x2(t) x3(t) x4(t) x5(t) x6(t)];
ode1 = subs(ode1,  [th, diff(th, t), diff(th, t, t),...
                    f, diff(f, t), diff(f, t, t),...
                    r, diff(r, t), diff(r, t, t)],...
                    [x1, x2, xth2, x3, x4, xf2, x5, x6, xr2]);
ode2 = subs(ode2,  [th, diff(th, t), diff(th, t, t),...
                    f, diff(f, t), diff(f, t, t),...
                    r, diff(r, t), diff(r, t, t)],...
                    [x1, x2, xth2, x3, x4, xf2, x5, x6, xr2]);
ode3 = subs(ode3,  [th, diff(th, t), diff(th, t, t),...
                    f, diff(f, t), diff(f, t, t),...
                    r, diff(r, t), diff(r, t, t)],...
                    [x1, x2, xth2, x3, x4, xf2, x5, x6, xr2]);

f1 = x2;
f2 = rhs(isolate(ode1, xth2));
f3 = x4;
f4 = rhs(isolate(ode2, xf2));
f5 = x6;
f6 = rhs(isolate(ode3, xr2));
f = [f1; f2; f3; f4; f5; f6];

%% linearization
syms x_1 x_2 x_3 x_4 x_5 x_6
x_ = [x_1 x_2 x_3 x_4 x_5 x_6];
f = expand(subs(f, x, x_));

J.symbolic = jacobian(f, x_);

syms x_1_0 x_2_0 x_3_0 x_4_0 x_5_0 x_6_0
x_0 = [x_1_0 x_2_0 x_3_0 x_4_0 x_5_0 x_6_0];

J.algebraic = simplify(subs(J.symbolic, x_, x_0));

% The solution requires a point where linearization occurs
point = [0, 0, 0, 0, 100, 1];
J.algebraic = simplify(subs(J.algebraic, x_0, point));
% J.algebraic

% linsys = J.algebraic * transpose(x - x_0);
linsys = J.algebraic * transpose(x - point);

eqns = subs(linsys, x, x_) == 0;

% Ax = b, X' = Ax - b (eqns = A*x_.' == 0)
[A, b] = equationsToMatrix(eqns, x_);

%% homogeneous system solution
[nA, mA] = size(A);
[nb, mb] = size(b);

if nA == nb
    n = nA;
else
    disp("Something wrong with converted matrix size");
    return
end

I = eye(n);
A = double(A);

% R is the matrix of right eigenvectors (A*R = R*D)
[R, D, L] = eig(A);

% d is the vector of eigenvalues
d = eig(A);

% u is the vector of unique eigenvalues
u = unique(d);

% X is a set of solutions
X = sym('X', n);

% fund is the matrix of fundamental system of solutions
fund = sym('fund', n);
% C is the vector of arbitrary constant numbers
C = sym('C', [1 n]);
% W is the Wronskian
W = sym('X', n);
% conj_arr is the array of complex conjugate
conj_arr = [];

for i = 1:length(u)
    % m is the algebraic multiplicity
    m = sum(d == d(i));
    % y is the geometric multiplicity
    y = n - rank(A - d(i)*I);
    if m == y && m == 1
        if isreal(d(i)) % if eigenvalue is real
            X(:,i) = exp(d(i)*t)*R(:,i);
        else % if eigenvalue is complex
            if ~ismember(d(i), conj_arr)
                % index of conjugate d(i)
                conj_i = find(d == conj(d(i)));
                % Euler's formula
                re_d = real(d(i));
                im_d = imag(d(i));
                c = cos(im_d*t);
                s = sin(im_d*t);
                re_v = real(R(:,i));
                im_v = imag(R(:,i));
                X(:,i) = exp(re_d*t)*(re_v*c - im_v*s);
                X(:,conj_i) = exp(re_d*t)*(im_v*c + re_v*s);
                conj_arr = [conj_arr, conj(d(i))];
            end
        end
    else
        disp("Something wrong with multiplicity!") % TODO (m == y && m > 1) & (y < m)
        return
    end
end

%% general homogeneous system solution

hss = X*C.';

% verifying result
%{
dhss = diff(hss, t);
hss_sys = dhss - A*hss;
disp("Verification hss")
v = vpa(simplify(hss_sys), 5);
pretty(v)
return
%}

%% particular case of the method of undetermined coefficients
A0 = sym('a0', [n 1]);
X1 = A0;
dX1 = diff(X1, t);
sys = dX1 == A*X1 - b;
sol = solve(sys);
sol = struct2cell(sol);
sol = [sol{1:6}];
sol = sol.';

%% general non-homogeneous system solution
nhss = hss + sol;

% verifying result
%{
dnhss = diff(nhss, t);
nhss_sys = dnhss - A*nhss + b;
disp("Verification nhss")
v = vpa(simplify(nhss_sys), 5);
pretty(v)
return
%}

%% search for unknown constants

sys = subs(nhss, t, 0) == point.';
pretty(vpa(sys, 3))
sol = solve(sys, C);
sol = struct2cell(sol);
sol = [sol{1:6}];
sol = sol.'

