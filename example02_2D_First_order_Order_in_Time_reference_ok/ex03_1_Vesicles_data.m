function pde = ex03_1_Vesicles_data(para)
if nargin == 0
    epsilon = 1;
    M = 1;
    C0 = 0;
    beta_m = 1;
    M1  = 0;
    M2 = 0;
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
    if ~isfield(para,'name') || isempty(para.name)
        name = 'ex03_1_Vesicles_data';
    else
        name = para.name;
    end
end

pde = struct('epsilon',epsilon, ...
    'M',M, ...
    'C0',C0, ...
    'beta_m', beta_m, ... 
    'M1',M1, ...
    'M2',M2, ...
    'init',@init, ...
    'name',name);

    function z = init(x,y)
        z = exact(x,y,0);
    end

    function z = exact(x,y,t)
%         z = exp(-t).^2.*sin(2.*x).*cos(2*y)./4;
        z = (sin(2*x).*sin(2*y)/4+0.48).*(1-sin(t).^2/2);
    end
end