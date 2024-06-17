function pde = ex03_Vesicles_data(para)
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
        name = 'ex03_Vesicles_data';
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
        z = (sin(2*x).*sin(2*y)/4+0.48).*(1-sin(t).^2/2);
    end

    function z = exact_t(x,y,t)
        z = -cos(t).*sin(t).*((sin(2.*x).*sin(2.*y))./4 + 12./25);
    end

    function z = rhs(x,y,t,beta)
        z = epsilon.*((4.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1) + 3.*cos(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3 - (9.*sin(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./4 - (3.*sin(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./20000 + (9.*cos(2.*x).^2.*sin(2.*x).*sin(2.*y).^3.*(sin(t).^2 - 2).^3)./8)./epsilon.^2 + (4.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1) + (3.*cos(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./2 + (3.*cos(2.*y).^2.*sin(2.*x).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./2 - (3.*sin(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./4 - (3.*sin(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./20000 - (3.*cos(2.*x).^2.*cos(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./2 + (3.*cos(2.*x).^2.*sin(2.*x).*sin(2.*y).^3.*(sin(t).^2 - 2).^3)./16 + (3.*cos(2.*y).^2.*sin(2.*x).^3.*sin(2.*y).*(sin(t).^2 - 2).^3)./16 - (3.*cos(2.*x).^2.*cos(2.*y).^2.*sin(2.*x).*sin(2.*y).*(sin(t).^2 - 2).^3)./4)./epsilon.^2 - 64.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)) + epsilon.*((4.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1) + 3.*cos(2.*y).^2.*sin(2.*x).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3 - (9.*sin(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./4 - (3.*sin(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./20000 + (9.*cos(2.*y).^2.*sin(2.*x).^3.*sin(2.*y).*(sin(t).^2 - 2).^3)./8)./epsilon.^2 + (4.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1) + (3.*cos(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./2 + (3.*cos(2.*y).^2.*sin(2.*x).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./2 - (3.*sin(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./4 - (3.*sin(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./20000 - (3.*cos(2.*x).^2.*cos(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./2 + (3.*cos(2.*x).^2.*sin(2.*x).*sin(2.*y).^3.*(sin(t).^2 - 2).^3)./16 + (3.*cos(2.*y).^2.*sin(2.*x).^3.*sin(2.*y).*(sin(t).^2 - 2).^3)./16 - (3.*cos(2.*x).^2.*cos(2.*y).^2.*sin(2.*x).*sin(2.*y).*(sin(t).^2 - 2).^3)./4)./epsilon.^2 - 64.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)) + (((3.*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^2)./40000 - 1).*((sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1) + (3.*cos(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./16 - (3.*sin(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./80000)./epsilon.^2 - 8.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)))./epsilon + (((3.*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^2)./40000 - 1).*((sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1) + (3.*cos(2.*y).^2.*sin(2.*x).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./16 - (3.*sin(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./80000)./epsilon.^2 - 8.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)))./epsilon - (cos(t).*sin(t).*((sin(2.*x).*sin(2.*y))./4 + 12./25))./M + M2.*epsilon.*((sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1) + (3.*cos(2.*x).^2.*sin(2.*y).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./16 - (3.*sin(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./80000)./epsilon.^2 - 8.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)).*((pi.^2.*(34978104032.*sin(t).^2 - 1040618024.*sin(t).^4 - 3851953992.*sin(t).^6 + 481494249.*sin(t).^8 + 59717987984))./(102400000000.*epsilon) - beta + (epsilon.*pi.^2.*(sin(t).^2 - 2).^2)./16) + M2.*epsilon.*((sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1) + (3.*cos(2.*y).^2.*sin(2.*x).^2.*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^3)./16 - (3.*sin(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./80000)./epsilon.^2 - 8.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)).*((pi.^2.*(34978104032.*sin(t).^2 - 1040618024.*sin(t).^4 - 3851953992.*sin(t).^6 + 481494249.*sin(t).^8 + 59717987984))./(102400000000.*epsilon) - beta + (epsilon.*pi.^2.*(sin(t).^2 - 2).^2)./16) + (3.*cos(2.*x).^2.*sin(2.*y).^2.*(sin(t).^2 - 2).^2.*((((25.*sin(2.*x).*sin(2.*y) + 48).^3.*(sin(t).^2 - 2).^3)./8000000 - (sin(t).^2./2 - 1).*((sin(2.*x).*sin(2.*y))./4 + 12./25))./epsilon.^2 + 2.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)))./(8.*epsilon) + (3.*cos(2.*y).^2.*sin(2.*x).^2.*(sin(t).^2 - 2).^2.*((((25.*sin(2.*x).*sin(2.*y) + 48).^3.*(sin(t).^2 - 2).^3)./8000000 - (sin(t).^2./2 - 1).*((sin(2.*x).*sin(2.*y))./4 + 12./25))./epsilon.^2 + 2.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)))./(8.*epsilon) - (3.*cos(2.*x).*sin(2.*y).*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(((cos(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1))./2 - (3.*cos(2.*x).*sin(2.*y).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./160000)./epsilon.^2 - 4.*cos(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)).*(sin(t).^2 - 2).^2)./(2.*epsilon) - (3.*cos(2.*y).*sin(2.*x).*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(((cos(2.*y).*sin(2.*x).*(sin(t).^2./2 - 1))./2 - (3.*cos(2.*y).*sin(2.*x).*(25.*sin(2.*x).*sin(2.*y) + 48).^2.*(sin(t).^2 - 2).^3)./160000)./epsilon.^2 - 4.*cos(2.*y).*sin(2.*x).*(sin(t).^2./2 - 1)).*(sin(t).^2 - 2).^2)./(2.*epsilon) - (3.*sin(2.*x).*sin(2.*y).*((sin(2.*x).*sin(2.*y))./4 + 12./25).*(sin(t).^2 - 2).^2.*((((25.*sin(2.*x).*sin(2.*y) + 48).^3.*(sin(t).^2 - 2).^3)./8000000 - (sin(t).^2./2 - 1).*((sin(2.*x).*sin(2.*y))./4 + 12./25))./epsilon.^2 + 2.*sin(2.*x).*sin(2.*y).*(sin(t).^2./2 - 1)))./epsilon;
    end
end