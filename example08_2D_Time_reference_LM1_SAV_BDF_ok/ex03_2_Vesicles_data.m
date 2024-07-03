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
    'name',name);

    function z = init(x,y)
        rl =[0.28*pi,0.28*pi];
        xl = [pi,pi];
        yl = [1.29*pi, -0.29*pi+pi];
        
        z =  tanh((rl(1)-sqrt((x-xl(1)).^2+(y-yl(1)).^2))./(sqrt(2)*epsilon)) ...
            +tanh((rl(2)-sqrt((x-xl(2)).^2+(y-yl(2)).^2))./(sqrt(2)*epsilon))+1;
    end
end