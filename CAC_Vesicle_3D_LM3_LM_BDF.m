function [phi,eta,lambda] = CAC_Vesicle_3D_LM3_LM_BDF(pde,domain,Nx,Ny,Nz,time,option)
% Solve 3D Vesicle model
% Qi Li
% 07/05/2024
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
    figname_mass = [pde.name,'_S1_',num2str(pde.S1),'_dt_',num2str(time.dt),'_mass.txt'];
    figname_energy = [pde.name,'_S1_',num2str(pde.S1),'_dt_',num2str(time.dt),'_energy.txt'];     
    out1 = fopen(figname_mass,'a');
    out2 = fopen(figname_energy,'a');
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

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

%% initialization phi 1
time1 = time; time1.T = dt;
[phi1,~,~] = CAC_Vesicle_3D_LM3_LM_1st(pde,domain,Nx,Ny,Nz,time1,option);

t = t+dt;

for nt = 2:nplot
    t = t+dt;
     
    phi_star = 2*phi1-phi0;  
    g1 = (4*phi1-phi0)/(2*dt*M)...
           + S1./epsilon.^3.*(phi_star - 1./(Lx.*Ly.*Lz).*fun_inner(phi_star,1)) ...
           - S2./epsilon.*(lap_diff(phi_star) - 1./(Lx.*Ly.*Lz).*fun_inner(lap_diff(phi_star),1)) ...
           + S3*epsilon.*(lap_diff(lap_diff(phi_star)) - 1./(Lx.*Ly.*Lz).*fun_inner(lap_diff(lap_diff(phi_star)),1));
    
    if isfield(pde,'rhs') && isfield(pde,'exact')
%         ephi   = pde.exact(xx,yy,t);
%         ephi_t = pde.exact_t(xx,yy,t);
%         emu   = epsilon*lap_diff(lap_diff(ephi)) + delta_Q(ephi);
%         tmp  = ephi_t./M - lap_diff(emu);
        rhs = pde.rhs(xx,yy,zz,t);
    else
        rhs = 0;
    end
    % step 1
    g1 = g1 + rhs;
    g2 = -(delta_Q(phi_star)-fun_inner(1,delta_Q(phi_star))./(Lx*Ly*Lz));
    g3 =   delta_B(phi_star)-fun_inner(1,delta_B(phi_star))./(Lx*Ly*Lz);

    phi1_new = inv_A(g1);
    phi2_new = inv_A(g2);
    phi3_new = inv_A(g3);

    % Step 3
    lambda = fun_inner(delta_B(phi_star),4*phi1-phi0-3*(phi1_new+phi2_new)) ./ fun_inner(delta_B(phi_star),3*phi3_new);
    lambda = fsolve(@(lambda)non_fun_lambda(lambda,phi1_new,phi2_new,phi3_new,phi00),lambda,opts);

    eta_initial = 1; 
    eta = fsolve(@(eta)non_fun_eta(eta,lambda,phi1_new,phi2_new,phi3_new,phi1,phi0),eta_initial,...
                   opts);
    
    % Step 4
    phi = phi1_new + eta.*phi2_new + lambda.*phi3_new; 
    
    %% update phi0
    phi0 = phi1; 
    phi1 = phi; 
    
    if 1 == option.energyflag
        calculate_energy(out1,out2,t,phi1,phi0)
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('lambda=%.4e,eta=%.4e,epsilon=%.3f,t=%.5f/%.4f, dt=%.2e, Nx=%d, Ny=%d, Nz=%d, timeElapsed=%f\n',lambda,eta,epsilon,t,T,dt,Nx,Ny,Nz,timeElapsed);
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

function result = inv_A(phi)
global dt k2 M S1 S2 S3 epsilon
    L1 = epsilon.*k2.^2;
    L2 = S1/epsilon.^3 - S2/epsilon*k2 + S3.*epsilon.*k2.^2;
    phihat = fftn(phi);
    r      = phihat./(3/2/dt/M + L1 + L2);
    r(1,1,1) = phihat(1,1,1)./(3/2/dt/M + L1(1,1,1) - L1(1,1,1)+ L2(1,1,1) - L2(1,1,1));
    result = real(ifftn(r));
end

function result = fun_Q(phi)
global epsilon
    result = 6./epsilon.^2.*phi.^2.*grad_square(phi) ...
             - 2/epsilon.^2.*grad_square(phi) ...
             + 1/epsilon.^4.*(f(phi)).^2;
    result = result*epsilon/2;
end

function result = delta_Q(phi)
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

function [] = calculate_energy(out1,out2,t,phi1,phi0)
global epsilon S1 S2 S3

energy_linear = fun_inner(1,1./2.*epsilon.*lap_diff(phi1).^2);
energy_nonlinear = fun_inner(1,fun_Q(phi1));
energy_original = energy_linear + energy_nonlinear;

phi_star = 2*phi1-phi0;

energy_S = S1./2./epsilon.^3.*(phi1 - phi0).^2 ...
         + S2./2./epsilon.*grad_square(phi1 - phi0) ...
         + S3.*epsilon./2.*lap_diff(phi1 - phi0).^2;
energy_S = fun_inner(1,energy_S);

energy_linear_bdf = fun_inner(1,epsilon/4.*(lap_diff(phi1).^2+lap_diff(phi_star).^2));
energy_phi1 = fun_inner(1,fun_Q(phi1));
energy_phi0 = fun_inner(1,fun_Q(phi0));
energy_original_discrete = energy_linear_bdf + energy_S + (3*energy_phi1 - energy_phi0)/2;

mass    = fun_inner(1,phi1);
surface = B(phi1);

fprintf(out1,'%14.6e  %.10f %.10f \n',t,mass,surface);
fprintf(out2,'%14.6e  %f  %f\n',t,energy_original,energy_original_discrete);
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

function r = non_fun_lambda(lambda,phi_new1,phi_new2,phi_new3,phi0)
    left = B(phi_new1+phi_new2+lambda.*phi_new3);
    right = B(phi0);
    r = left - right;
end

function r = non_fun_eta(eta,lambda,phi_new1,phi_new2,phi_new3,phi1,phi0)
    phi_star = 2*phi1 - phi0 ;
    left = fun_inner(1,3.*fun_Q(phi_new1+eta.*phi_new2+lambda.*phi_new3)-4.*fun_Q(phi1)+fun_Q(phi0));
    right = eta.*fun_inner(delta_Q(phi_star),3.*(phi_new1+eta.*phi_new2+lambda.*phi_new3)-4.*phi1+phi0)...
             -lambda.*fun_inner(delta_B(phi_star),3.*(phi_new1+eta.*phi_new2+lambda.*phi_new3)-4.*phi1+phi0);
    r = left - right;
end

