function pde = ex03_2_Vesicles_data(para)
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
        name = 'ex03_2_Vesicles_data';
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
    'exact',@exact, ...
    'exact_t',@exact_t, ...
    'rhs',@rhs, ...
    'name',name);

    function z = init(x,y)
        z = exact(x,y,0);
    end

    function z = exact(x,y,t)
        z = sin(x).*sin(y).*sin(t);
    end

    function z = exact_t(x,y,t)
        z = cos(t).*sin(x).*sin(y);
    end

    function z = rhs(x,y,t,beta)
        z = (sin(x).*sin(y).*(64.*M.*sin(t) + 32.*epsilon.^3.*cos(t) - 256.*M.*epsilon.^2.*sin(t) + 256.*M.*epsilon.^4.*sin(t) + 768.*M.*cos(x).^2.*sin(t).^3.*sin(y).^2 + 768.*M.*cos(y).^2.*sin(t).^3.*sin(x).^2 - 768.*M.*sin(t).^3.*sin(x).^2.*sin(y).^2 + 960.*M.*sin(t).^5.*sin(x).^4.*sin(y).^4 + 64.*M.*M2.*beta.*epsilon.^2.*sin(t) - 128.*M.*M2.*beta.*epsilon.^4.*sin(t) + 2304.*M.*epsilon.^2.*cos(x).^2.*cos(y).^2.*sin(t).^3 - 4224.*M.*epsilon.^2.*cos(x).^2.*sin(t).^3.*sin(y).^2 - 4224.*M.*epsilon.^2.*cos(y).^2.*sin(t).^3.*sin(x).^2 + 3072.*M.*epsilon.^2.*sin(t).^3.*sin(x).^2.*sin(y).^2 + 32.*M.*M2.*epsilon.*pi.^2.*sin(t).^3 + 128.*M.*M2.*epsilon.^3.*pi.^2.*sin(t) - 9.*M.*M2.*epsilon.*pi.^2.*sin(t).^5 - 128.*M.*M2.*epsilon.^3.*pi.^2.*sin(t).^3 + 18.*M.*M2.*epsilon.^3.*pi.^2.*sin(t).^5 + 128.*M.*M2.*epsilon.^5.*pi.^2.*sin(t).^3 - 1920.*M.*cos(x).^2.*sin(t).^5.*sin(x).^2.*sin(y).^4 - 1920.*M.*cos(y).^2.*sin(t).^5.*sin(x).^4.*sin(y).^2 - 64.*M.*M2.*epsilon.*pi.^2.*sin(t) - 192.*M.*M2.*epsilon.^3.*pi.^2.*cos(x).^2.*sin(t).^5.*sin(y).^2 - 192.*M.*M2.*epsilon.^3.*pi.^2.*cos(y).^2.*sin(t).^5.*sin(x).^2 + 192.*M.*M2.*epsilon.^3.*pi.^2.*sin(t).^5.*sin(x).^2.*sin(y).^2 - 192.*M.*M2.*epsilon.*pi.^2.*cos(x).^2.*sin(t).^3.*sin(y).^2 - 192.*M.*M2.*epsilon.*pi.^2.*cos(y).^2.*sin(t).^3.*sin(x).^2 + 96.*M.*M2.*epsilon.*pi.^2.*cos(x).^2.*sin(t).^5.*sin(y).^2 + 96.*M.*M2.*epsilon.*pi.^2.*cos(y).^2.*sin(t).^5.*sin(x).^2 - 27.*M.*M2.*epsilon.*pi.^2.*cos(x).^2.*sin(t).^7.*sin(y).^2 - 27.*M.*M2.*epsilon.*pi.^2.*cos(y).^2.*sin(t).^7.*sin(x).^2 + 192.*M.*M2.*epsilon.*pi.^2.*sin(t).^3.*sin(x).^2.*sin(y).^2 - 96.*M.*M2.*epsilon.*pi.^2.*sin(t).^5.*sin(x).^2.*sin(y).^2 + 27.*M.*M2.*epsilon.*pi.^2.*sin(t).^7.*sin(x).^2.*sin(y).^2 + 192.*M.*M2.*beta.*epsilon.^2.*cos(x).^2.*sin(t).^3.*sin(y).^2 + 192.*M.*M2.*beta.*epsilon.^2.*cos(y).^2.*sin(t).^3.*sin(x).^2 - 192.*M.*M2.*beta.*epsilon.^2.*sin(t).^3.*sin(x).^2.*sin(y).^2))./(32.*M.*epsilon.^3);
    end
end