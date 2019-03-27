function NPV  = NPV(x, varargin)
    global NFE
    if isempty(NFE)
        NFE=0;
    end
    NFE=NFE+1;

    % Objective function should be prepared to receive a set of points

    % From Optimization Algorithm we get x
    % x is between 0 and 1

    %... here we have a pattern , and the center of pattern is goal point
    %... this point cannot be at corner 

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
    ECLDATA.CompLoc(:, 3) = 1;
    ECLDATA.CompLoc(:, 4) = 3;

    ECLDATA=Eclipse(ECLDATA, varargin{3});

    Qopcum=ECLDATA.Qop;    % This is FGPT
    Qwpcum=ECLDATA.Qwp;    % This FWPT
    Qwicum=ECLDATA.Qwi;   % This is cumulative  FWIT
    Qgpcum=ECLDATA.Qgp;
    Qopcum=Qopcum';
    Qwpcum=Qwpcum';
    Qwicum=Qwicum';
    Qgpcum=Qgpcum';

    Qop=Qopcum;Qwp=Qwpcum; Qwi=Qwicum; Qgp=Qgpcum;
    %now we want to obtain rate from cumulative prod
    for i=2:numel(Qopcum)
        Qop(i)=Qop(i)-Qopcum(i-1);
        Qwp(i)=Qwp(i)-Qwpcum(i-1);
        Qwi(i)=Qwi(i)-Qwicum(i-1);                 % Water injected STB at every time step
        Qgp(i)=Qgp(i)-Qgpcum(i-1);
    end
    
    t=ECLDATA.Time;
    t=(t');


    % Simplified form of NPV
    Nwell=numel(x)/2; % + varargin{2}; % number of wells
    % Economic Parameters :
    Ccomp=5e6;                 %Each Well compeletion cost ($)
    Cdrill=1000/0.3048;        % Drilling Cost ($/m)
    Po=80*6.287;               % Oil price ($/SM3)
    Pg=240;                    % gas price ($/MSM3)
    Pwi=8*6.287;               % Water injected cost ($/SM3)
    Pwp=10*6.287;              % Water Produced cost ($/SM3)
    r=0.1;                     % Discount Rate
    Lwmain=2717.9016;          % Length of the main bore (m)

    % Capital Expenditure (total cost of Drill& Complete well) ($)
    Ccapex=Nwell*(Ccomp+(Lwmain).*Cdrill);

    Qgp=Qgp/(10^3); %Convert to MSM3
    Rt=Po.*Qop+Pg.*Qgp;    % Renevue  at time t
    Et=Pwp.*Qwp+Pwi.*Qwi;
    CFt=Rt-Et;               % Cash Flow at time t

    %%This is Net Present Value: (Revenue - Cost)/discount - capex
    NPV=(sum(CFt./(1+r).^t)-Ccapex);

return


