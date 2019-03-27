function NPV  = NPV(x, varargin)

    global NFE
    % NFE=0;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1

    % Objective function should be prepared to receive a set of points
    % appropriate for rectangular serahc domain

    % From Optimization Algorithm we get x
    % x is between 0 and 1
    %... here we have a pattern , and the center of pattern is goal point
    %.. this point cannot be at corner 
    %
    nwi = numel(x)/2;   %number of wells to be optimized
    if(nwi > 1)
        for i=1:2
            if(mod(i,2) ~= 0)
                x(i:2:end) = min(floor(1+(varargin{1}).*x(i:2:end)), varargin{1});
            else
                x(i:2:end) = min(floor(1+(varargin{2}).*x(i:2:end)), varargin{2});
            end
        end
    else
        x(1) = min(floor(1+(varargin{1})*x(1)), varargin{1});
        x(2) = min(floor(1+(varargin{1})*x(2)), varargin{2});
    end
    
    for i=1:nwi
        ECLDATA.CompLoc(i,1:2) = x(2*i-1:2*i);
        ECLDATA.WHLoc(i,:) = x(2*i-1:2*i);
    end
    ECLDATA.CompLoc(:, 3:4)=1;

    ECLDATA=Eclipse(ECLDATA, varargin{3});

    Qopcum=ECLDATA.Qop;    % This is FGPT
    Qwpcum=ECLDATA.Qwp;    % This FWPT
    Qwicum=ECLDATA.Qwi;   % This is cumulative  FWIT
    Qopcum=Qopcum';
    Qwpcum=Qwpcum';
    Qwicum=Qwicum';

    Qop=Qopcum;Qwp=Qwpcum; Qwi=Qwicum;
    %now we want to obtain rate from cumulative prod
    for i=2:numel(Qopcum)
        Qop(i)=Qop(i)-Qopcum(i-1);
        Qwp(i)=Qwp(i)-Qwpcum(i-1);
        Qwi(i)=Qwi(i)-Qwicum(i-1);                 % Water injected STB at every time step
    end

    t=ECLDATA.Time;
    t=(t');


    % Simplified form of NPV
    Nwell=numel(x)/2; % + varargin{3}; % number of wells 
    % Economic Parameters :
    Ccomp=5e6;            % Each well compeletion cost ($)
    Cdrill=1000;          % Drilling Cost ($/ft)
    Po=80;                % Oil price ($/STB)
    Pg=240*0.3048^3;      % gas price ($/MScf)
    Pwi=8;                % Water injected cost ($/STB)
    Pwp=10;               % Water Produced cost ($/STB)
    r=0.1;                % Discount Rate
    Lwmain=8917;          % Length of the main bore (m)/0.3048=ft 

    % Capital Expenditure (total cost of Drill& Complete well) ($)
    Ccapex=Nwell*(Ccomp +(Lwmain).*Cdrill);

    Qgp=0;
    Rt=Po.*Qop+Pg.*Qgp;    % Renevue  at time t
    Et=Pwp.*Qwp+Pwi.*Qwi;
    CFt=Rt-Et;               % Cash Flow at time t

    %%This is Net Present Value: (Revenue - Cost)/discount - capex
    NPV=(sum(CFt./(1+r).^t)-Ccapex);

return


