%% Epitaxy
clear 
clc
syms x y epsilon t M alpha beta

% phi = (sin(2*x).*sin(2*y)/4+0.48).*(1-sin(t).^2/2);
phi = sin(x).*sin(y).*sin(t);
lap_phi = diff(phi,x,2) + diff(phi,y,2);
grad_square = diff(phi,x,1).^2 + diff(phi,y,1).^2;

F     = (phi.^2-1).^2/4;
f     = phi.^3 - phi;
f_der = 3*phi.^2 - 1;

omega = -lap_phi +1./epsilon.^2.*f;
lap_omega = diff(omega,x,2) + diff(omega,y,2);

mu = -epsilon.*lap_omega + 1./epsilon.*f_der.*omega;
lap_mu = diff(mu,x,2) + diff(mu,y,2);

f1 = simplify( diff(phi,t)./M - lap_mu);
f2 = simplify( diff(phi,t));

f1 = char(f1);
f1 = strrep(f1,'*','.*');
f1 = strrep(f1,'/','./');
f1 = strrep(f1,'^','.^')

f2 = char(f2);
f2 = strrep(f2,'*','.*');
f2 = strrep(f2,'/','./');
f2 = strrep(f2,'^','.^')


