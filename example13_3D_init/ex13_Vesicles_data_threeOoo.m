function pde = ex13_Vesicles_data_threeOoo(para)
if nargin == 0
    epsilon = 1;
    M = 1;
    C0 = 0;
    beta_m = 1;
    M1  = 0;
    M2 = 0;
    S1 = 0;
    S2 = 0;
    S3 = 0;
    name = '';
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'epsilon') || isempty(para.epsilon)
        epsilon = 1;
    else
        epsilon = para.epsilon;
    end
    if ~isfield(para,'M') || isempty(para.M)
        M = 1;
    else
        M = para.M;
    end
    if ~isfield(para,'C0') || isempty(para.C0)
        C0 = 0;
    else
        C0 = para.C0;
    end
    if ~isfield(para,'beta_m') || isempty(para.beta_m)
        beta_m = 1;
    else
        beta_m = para.beta_m;
    end
    if ~isfield(para,'M1') || isempty(para.M1)
        M1 = 0;
    else
        M1 = para.M1;
    end
    if ~isfield(para,'M2') || isempty(para.M2)
        M2 = 0;
    else
        M2 = para.M2;
    end
    if ~isfield(para,'S1') || isempty(para.S1)
        S1 = 0;
    else
        S1 = para.S1;
    end
    if ~isfield(para,'S2') || isempty(para.S2)
        S2 = 0;
    else
        S2 = para.S2;
    end
    if ~isfield(para,'S3') || isempty(para.S3)
        S3 = 0;
    else
        S3 = para.S3;
    end
    if ~isfield(para,'name') || isempty(para.name)
        name = 'ex13_Vesicles_data_threeOoo';
    else
        name = para.name;
    end
end

pde = struct( ...
    'epsilon',epsilon, ...
    'M',M, ...
    'C0',C0, ...
    'beta_m', beta_m, ... 
    'M1',M1, ...
    'M2',M2, ...
    'S1',S1, ...
    'S2',S2, ...
    'S3',S3, ...
    'init',@init, ...
    'exact',@exact, ...
    'name',name);

    function z = init(x,y,z)
        d = 0.152*pi;
        rl =[0.15*pi,0.15*pi,0.35*pi];
        xl = [pi,pi,pi];
        yl = [pi,pi,pi];
        zl = [pi+3.6*d,pi+1.2*d,pi-2.5*d];
        rx = 1;
        ry = 1;
        rz = 1;
        
        z =  tanh((rl(1)-sqrt((x-xl(1)).^2./rx.^2+(y-yl(1)).^2./ry.^2+(z-zl(1)).^2./rz.^2))./(sqrt(2)*epsilon)) ...
            +tanh((rl(2)-sqrt((x-xl(2)).^2./rx.^2+(y-yl(2)).^2./ry.^2+(z-zl(2)).^2./rz.^2))./(sqrt(2)*epsilon)) ...
            +tanh((rl(3)-sqrt((x-xl(3)).^2./rx.^2+(y-yl(3)).^2./ry.^2+(z-zl(3)).^2./rz.^2))./(sqrt(2)*epsilon)) +2;
    end
end