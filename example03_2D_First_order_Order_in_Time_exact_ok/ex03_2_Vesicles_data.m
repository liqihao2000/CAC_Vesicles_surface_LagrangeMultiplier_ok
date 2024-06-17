function pde = ex03_2_Vesicles_data(para)
if nargin == 0
    epsilon = 1;
    M = 1;
    C0 = 0;
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
    if ~isfield(para,'name') || isempty(para.name)
        name = 'ex03_2_Vesicles_data';
    else
        name = para.name;
    end
end

pde = struct('epsilon',epsilon, ...
    'M',M, ...
    'C0',C0, ...
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

    function z = rhs(x,y,t)
        z = (sin(x).*sin(y).*(8.*M.*sin(t) + epsilon.^3.*cos(t) + 88.*M.*epsilon.^2.*sin(t) + 8.*M.*epsilon.^4.*sin(t) - 36.*M.*cos(t).^2.*sin(t) + 30.*M.*cos(t).^4.*sin(t) - 72.*M.*cos(x).^2.*sin(t) + 90.*M.*cos(x).^4.*sin(t) - 72.*M.*cos(y).^2.*sin(t) + 90.*M.*cos(y).^4.*sin(t) + 192.*M.*cos(t).^2.*cos(x).^2.*sin(t) - 180.*M.*cos(t).^2.*cos(x).^4.*sin(t) - 120.*M.*cos(t).^4.*cos(x).^2.*sin(t) + 90.*M.*cos(t).^4.*cos(x).^4.*sin(t) + 192.*M.*cos(t).^2.*cos(y).^2.*sin(t) - 180.*M.*cos(t).^2.*cos(y).^4.*sin(t) - 120.*M.*cos(t).^4.*cos(y).^2.*sin(t) + 90.*M.*cos(t).^4.*cos(y).^4.*sin(t) + 288.*M.*cos(x).^2.*cos(y).^2.*sin(t) - 240.*M.*cos(x).^2.*cos(y).^4.*sin(t) - 240.*M.*cos(x).^4.*cos(y).^2.*sin(t) + 150.*M.*cos(x).^4.*cos(y).^4.*sin(t) - 96.*M.*epsilon.^2.*cos(t).^2.*sin(t) - 228.*M.*epsilon.^2.*cos(x).^2.*sin(t) - 228.*M.*epsilon.^2.*cos(y).^2.*sin(t) - 648.*M.*cos(t).^2.*cos(x).^2.*cos(y).^2.*sin(t) + 480.*M.*cos(t).^2.*cos(x).^2.*cos(y).^4.*sin(t) + 480.*M.*cos(t).^2.*cos(x).^4.*cos(y).^2.*sin(t) + 360.*M.*cos(t).^4.*cos(x).^2.*cos(y).^2.*sin(t) - 300.*M.*cos(t).^2.*cos(x).^4.*cos(y).^4.*sin(t) - 240.*M.*cos(t).^4.*cos(x).^2.*cos(y).^4.*sin(t) - 240.*M.*cos(t).^4.*cos(x).^4.*cos(y).^2.*sin(t) + 150.*M.*cos(t).^4.*cos(x).^4.*cos(y).^4.*sin(t) + 228.*M.*epsilon.^2.*cos(t).^2.*cos(x).^2.*sin(t) + 228.*M.*epsilon.^2.*cos(t).^2.*cos(y).^2.*sin(t) + 432.*M.*epsilon.^2.*cos(x).^2.*cos(y).^2.*sin(t) - 432.*M.*epsilon.^2.*cos(t).^2.*cos(x).^2.*cos(y).^2.*sin(t)))./(M.*epsilon.^3);
    end
end