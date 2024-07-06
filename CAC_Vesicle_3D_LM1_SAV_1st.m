function [phi, u0] = CAC_Vesicle_3D_LM1_SAV_1st(pde,domain,Nx,Ny,Nz,time,option)
% Solve 3D Vesicle model
% Qi Li
% 07/04/2024
global dt kx ky kz kxx kyy kzz k2 k4 hx hy hz Lx Ly Lz ...
       epsilon M S1 S2 S3 C0

if ~exist('option','var'), option = []; end
if ~isfield(option,'tol')
    option.tol = 10^-14;   % default tol
end
if ~isfield(option,'tolit')
    option.tolit = 10^-8;   % default tolit
end
if ~isfield(option,'maxit')
    option.maxit = 2000;   % default maxit
end
if ~isfield(option,'plotflag')
    option.plotflag = 0;   
end
if ~isfield(option,'saveflag')
    option.saveflag = 0;  
end
if ~isfield(option,'savefinal')
    option.savefinal = 0;  
end
if ~isfield(option,'printflag')
    option.printflag = 0;   
end
if ~isfield(option,'vtkflag')
    option.printflag = 0;   
end
if ~isfield(option,'energyflag')
    option.energyflag = 0;   
end
if 1 == option.energyflag
    figname_mass = [pde.name,num2str(time.dt),'_mass.txt'];
    figname_energy = [pde.name,num2str(time.dt),'_energy.txt'];      
    out1 = fopen(figname_mass,'w');
    out2 = fopen(figname_energy,'w');
end

%% options for fsolve
opts = optimoptions('fsolve','Display','off',...
                       'StepTolerance', option.tol, ...
                       'FunctionTolerance',option.tolit,...
                       'MaxIterations',option.maxit);


%%
T  = time.T;
t  = time.t0;
dt = time.dt;
tsave = time.tsave;

dir_fig  = [pde.name '/fig'];
dir_data = [pde.name '/data'];

epsilon = pde.epsilon;
M       = pde.M;
S1      = pde.S1;
S2      = pde.S2;
S3      = pde.S3;
C0      = pde.C0;

Lx = domain.xright - domain.xleft;
Ly = domain.yright - domain.yleft;
Lz = domain.zright - domain.zleft;
hx = Lx/Nx;
hy = Ly/Ny;
hz = Lz/Nz;
x  = domain.xleft  + hx*(0:Nx-1);
y  = domain.yleft  + hy*(0:Ny-1);
z  = domain.zleft  + hz*(0:Nz-1);

% Fourier spectral 
k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);
k_z = 1i*[0:Nz/2 -Nz/2+1:-1]*(2*pi/Lz);
% [kx, ky, kz] = ndgrid(k_x,k_y,k_z);
k2x = k_x.^2;
k2y = k_y.^2;
k2z = k_z.^2;
[kxx, kyy, kzz] = ndgrid(k2x,k2y,k2z);
k2 = kxx + kyy + kzz;
k4 = k2.^2;
% Furthermore, it is important to make the highest frequency N/2 to zero 
% in odd derivatives due to the symmetry.
k_x1 = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]*(2*pi/Lx);
k_y1 = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]*(2*pi/Ly);
k_z1 = 1i*[0:Nz/2-1 0 -Nz/2+1:-1]*(2*pi/Lz);
[kx, ky, kz] = ndgrid(k_x1,k_y1,k_z1);


[xx,yy,zz] = ndgrid(x,y,z);
phi0 = pde.init(xx,yy,zz);
phi00 = phi0;
nfigure =1;

%% plot initial value
nfigure =1;
if 1 == option.saveflag
    if ~exist(dir_data,'dir')
        mkdir(dir_data);
    end
    ss1 = [dir_data '/phi_t=' num2str(t) '.txt'];
    writematrix(phi0,ss1,'Delimiter',' ');
    writematrix(xx,[dir_data '/X.txt'],'Delimiter',' ');
    writematrix(yy,[dir_data '/Y.txt'],'Delimiter',' ');
    writematrix(zz,[dir_data '/Y.txt'],'Delimiter',' ');
end
if 1 == option.plotflag
    if 1 == option.saveflag
        showsolution_2D(nfigure,xx,yy,phi0,t,dir_fig);
    else
        showsolution_2D(nfigure,xx,yy,phi0,t);
    end
end

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

u0 = fun_u_init(phi0);

% Initial energy
if 1 == option.energyflag
    calculate_energy1(out1,out2,hx,hy,hz,t,phi0,r0);
end

for nt = 1:nplot
    t = t+dt;
     
    phi_star = phi0;    
    
    % step 1
    H = fun_H(phi_star);
    g = u0 - 1/2*fun_inner(H,phi0);
    
    g1 = phi0/dt/M...
           + S1./epsilon.^3.*(phi_star - 1./(Lx.*Ly.*Lz).*fun_inner(phi_star,1)) ...
           - S2./epsilon.*(lap_diff(phi_star) - 1./(Lx.*Ly.*Lz).*fun_inner(lap_diff(phi_star),1)) ...
           + S3*epsilon.*(lap_diff(lap_diff(phi_star)) - 1./(Lx.*Ly.*Lz).*fun_inner(lap_diff(lap_diff(phi_star)),1))...
           - H.*g;
    if isfield(pde,'rhs') && isfield(pde,'exact')
%         ephi   = pde.exact(xx,yy,t);
%         ephi_t = pde.exact_t(xx,yy,t);
%         emu   = epsilon*lap_diff(lap_diff(ephi)) + fun_q(ephi);
%         tmp  = ephi_t./M - lap_diff(emu);
        tmp = pde.rhs(xx,yy,zz,t);
        g1 = g1 + tmp;
    end
    g2 = delta_B(phi_star)-fun_inner(1,delta_B(phi_star))./(Lx*Ly*Lz);

    psiA = inv_A(H);
    psiB = inv_A(g1);  
    psiC = inv_A(g2); 

    b1 = fun_inner(H, psiB);
    b2 = fun_inner(H, psiC);
    
    mu = 0.5*fun_inner(H,psiA);
    
    X = b1./(mu + 1);
    Y = b2./(mu + 1);
   
    % Step 3
    phi1_new = -1/2*X.*psiA + psiB;
    phi2_new = -1/2*Y.*psiA + psiC;
    
    u1_new = g + 1/2*X;
    u2_new = 1/2*Y;

    % Step 3
    lambda = fun_inner(delta_B(phi_star),phi0-phi1_new) ./ fun_inner(delta_B(phi_star),phi2_new);

    lambda = fsolve(@(lambda)non_fun(lambda,phi1_new,phi2_new,phi00),lambda,...
                   opts);

    phi = phi1_new + lambda*phi2_new; 
    u   = u1_new   + lambda*u2_new; 
    
    %% update phi0
    phi0 = phi;  
    u0 = u;
    
    if 1 == option.energyflag
        calculate_energy1(out1,out2,hx,hy,hz,t,phi0,r0);
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('lambda=%.4e,epsilon=%.3f,t=%.5f/%.4f, dt=%.2e, Nx=%d, Ny=%d, Nz=%d, timeElapsed=%f\n',lambda,epsilon,t,T,dt,Nx,Ny,Nz,timeElapsed);
        end
        
        if 1 == option.saveflag
            if ~exist(dir_data,'dir')
                mkdir(dir_data);
            end
            ss1 = [dir_data '/phi_t=' num2str(t) '.txt'];
            writematrix(phi0,ss1,'Delimiter',' ');
        end
        
        nfigure = nfigure +1;
        if 1 == option.plotflag
            if 1 == option.vtkflag
                write_vtk_grid_values(dir_data,x,y,nt,phi0);
            end
            if 1 == option.saveflag
                showsolution_2D(nfigure,xx,yy,phi,t,dir_fig);
            else
                showsolution_2D(nfigure,xx,yy,phi,t);
            end
        end
    end
    
end

if 1 == option.savefinal
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),'S1=',num2str(pde.S1),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt)];    
    filename=[name '.mat'];
    save(filename,'epsilon','xx','yy','zz','hx','hy','hz','Nx','Ny','Nz','dt','T','phi','domain');
end

if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end

function result = fun_u_init(phi)
global C0
if fun_inner(1,fun_Q(phi)) + C0 <0
    disp("Root < 0");
    return;
end
result  = sqrt(fun_inner(1,fun_Q(phi)) + C0);
end

function result = fun_H(phi)
global C0
if fun_inner(1,fun_Q(phi)) + C0 <0
    disp("Root < 0");
    return;
end
result  = fun_q(phi)./sqrt(fun_inner(1,fun_Q(phi)) + C0);
end

function result = inv_A(phi)
global dt k2 M epsilon Lx Ly Lz S1 S2 S3
    L1 = epsilon.*k2.^2;
    L2 = S1/epsilon.^3 - S2/epsilon*k2 + S3*epsilon.*k2.^2;
    phihat = fftn(phi);
    r      = phihat./(1/dt/M + L1 + L2);
    r(1,1,1) = phihat(1,1,1)./(1/dt/M + L1(1,1,1) + L2(1,1,1) - 1/(Lx*Ly*Lz).*L1(1,1,1).*Lx.*Ly.*Lz - 1/(Lx*Ly*Lz).*L2(1,1,1).*Lx.*Ly.*Lz);
    result = real(ifftn(r));
end

function result = fun_Q(phi)
global epsilon
    result = 6./epsilon.^2.*phi.^2.*grad_square(phi) ...
             - 2/epsilon.^2.*grad_square(phi) ...
             + 1/epsilon.^4.*(f(phi)).^2;
    result = result*epsilon/2;
end

function result = fun_q(phi)
global epsilon
    div_term = diff_x(phi.^2.*diff_x(phi)) + diff_y(phi.^2.*diff_y(phi))+ diff_z(phi.^2.*diff_z(phi));
    result = 6./epsilon.*(phi.*grad_square(phi) - div_term) ...
             + 2/epsilon.*lap_diff(phi) ...
             + 1/epsilon.^3.*f(phi).*f_der(phi);
end

function r = fun_inner(f,g)
global hx hy hz
    r1 = fftn(f.*g);
    r = r1(1,1,1)*hx*hy*hz;
end

function [] = calculate_energy1(out1,out2,hx,hy,hz,t,phi,r)
global C0
energy_linear = fftn(energyoperatorL(phi));
energy_nonlinear = fftn(F(phi));

energy_original = energy_linear(1,1,1)*hx*hy*hz + energy_nonlinear(1,1,1)*hx*hy*hz;

energy_modified = energy_linear(1,1,1)*hx*hy*hz + r.^2 - C0;

mass    = fftn(phi);
mass    = mass(1,1,1)*hx*hy*hz;

fprintf(out1,'%14.6e  %.8f \n',t,mass);
fprintf(out2,'%14.6e  %f  %f\n',t,energy_original,energy_modified);
end

function lap=lap_diff(phi)
global k2
    lap=real(ifftn((k2.*fftn(phi))));
end

function lap=diff_x(phi)
global kx
    lap=real(ifftn((kx.*fftn(phi))));
end

function lap=diff_y(phi)
global ky
    lap=real(ifftn((ky.*fftn(phi))));
end

function lap=diff_z(phi)
global kz
    lap=real(ifftn((kz.*fftn(phi))));
end

function result = grad_square(phi)
    result = diff_x(phi).^2 + diff_y(phi).^2 + diff_z(phi).^2;
end

function result = f_der(phi)
    result = 3.*phi.^2-1;
end

function result = f(phi)
    result = phi.^3 - phi;
end

function result = F(phi)
    result = 1/4*(phi.^2-1).^2;
end

function result = B(phi)
global hx hy hz epsilon
    r = fftn(epsilon./2.*grad_square(phi)+1/epsilon.*F(phi));
    result = r(1,1,1)*hx*hy*hz;
end

function result = delta_B(phi)
global epsilon
    result = -epsilon*lap_diff(phi) + 1/epsilon.*f(phi);
end

function r = non_fun(alpha,phi_new1,phi_new2,phi0)
    left = B(phi_new1+alpha.*phi_new2);
    right = B(phi0);
    r = left - right;
end


