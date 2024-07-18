%% Epitaxy
clear 
clc
syms x y epsilon t M Lx Ly alpha beta M1 M2 

% phi = (sin(2*x).*sin(2*y)/4+0.48).*(1-sin(t).^2/2);
% phi = (sin(2*x).*sin(2*y)/4).*(1-sin(t).^2/2);
phi = sin(x).*sin(y).*sin(t);
lap_phi = diff(phi,x,2) + diff(phi,y,2);
grad_square = diff(phi,x,1).^2 + diff(phi,y,1).^2;

F     = (phi.^2-1).^2/4;
f     = phi.^3 - phi;
f_der = 3*phi.^2 - 1;

omega = -lap_phi +1./epsilon.^2.*f;
lap_omega = diff(omega,x,2) + diff(omega,y,2);

% A_phi = int(int(phi+1,x,0,2*pi),y,0,2*pi);
B_phi = int(int(epsilon./2.*grad_square+1./epsilon.*F,x,0,2*pi),y,0,2*pi);

mu = -epsilon.*lap_omega + 1./epsilon.*f_der.*omega ...
     + M2.*(B_phi-beta).*epsilon.*omega;

fprintf("mubar\n")
mubar = mu - 1./(Lx*Ly)*int(int(mu,x,0,Lx),y,0,Ly);
f1 = simplify( diff(phi,t)./M + mubar)./epsilon;
f2 = simplify( diff(phi,t));

f1 = char(f1);
f1 = strrep(f1,'*','.*');
f1 = strrep(f1,'/','./');
f1 = strrep(f1,'^','.^')

f2 = char(f2);
f2 = strrep(f2,'*','.*');
f2 = strrep(f2,'/','./');
f2 = strrep(f2,'^','.^')

