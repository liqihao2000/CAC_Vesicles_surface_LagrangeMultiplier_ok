function pde = ex12_Vesicles_data_eightO(para)
if nargin == 0
    epsilon = 1;
    M = 1;
    C0 = 0;
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
        name = 'ex12_Vesicles_data_eightO';
    else
        name = para.name;
    end
end

pde = struct('epsilon',epsilon, ...
    'M',M, ...
    'C0',C0, ...
    'S1',S1, ...
    'S2',S2, ...
    'S3',S3, ...
    'init',@init, ...
    'exact',@exact, ...
    'name',name);

    function z = init(x,y)
        d = 0.31*pi;
        rl =[0.2*pi,0.2*pi,0.2*pi,0.2*pi,0.2*pi,0.2*pi,0.2*pi,0.2*pi];
        xl = [pi+0,pi+d,pi+0,pi-d,pi+d,  pi+2*d,pi-d,  pi-2*d];
        yl = [pi+d,pi+0,pi-d,pi+0,pi+2*d,pi-d,  pi-2*d,pi+d];
        rx = 1;
        ry = 1;
        
        z = -1;
        for k = 1:length(rl)
            z = z + tanh((rl(k)-sqrt((x-xl(k)).^2./rx.^2+(y-yl(k)).^2./ry.^2))./(sqrt(2)*epsilon)) + 1;
        end
    end
end