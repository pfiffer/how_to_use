function NPV  = NPV(x, varargin)

    global NFE
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    
    if(~isempty(varargin))
        SD = varargin{1};
        nx = numel(SD(:,1));
        ny = numel(SD(1,:));
    end
    
    % Objective function should be prepared to receive a set of points
    
    % From Optimization Algorithm we get x
    % x is between 0 and 1
    nwi = numel(x)/2;   % Injection well number
    if(nwi>1)
        x(1:2:end) = min(floor(1+nx.*x(1:2:end)), nx);
        x(2:2:end) = min(floor(1+ny.*x(2:2:end)), ny);
        k = SD(x(1:2:end), :);
    else
        x(1) = min(floor(1+(nx)*x(1)), nx);
        x(2) = min(floor(1+(ny)*x(2)), ny);
        % Selecting a row from Search Domain
        k = SD(x(1), :);
    end
    
    % As k has more than one rows and always has two columns, we should
    % change the code in order to math the dimension. In only PSO, x was a
    % row matrix and after changing the elements by the above loop, it
    % conserved it's shape
    for i=1:nwi
        ECLDATA.CompLoc(i,1:2) = k(i,1:2);
        ECLDATA.WHLoc(i,:) = k(i,1:2);
    end
    
    
    ECLDATA.CompLoc(:, 3) = 1;   %Completeion layers
    ECLDATA.CompLoc(:, 4) = 3;


    ECLDATA=Eclipse(ECLDATA, varargin{2});

    Qopcum=ECLDATA.Qop;    % This is FGPT
    Qwpcum=ECLDATA.Qwp;    % This FWPT
    Qwicum=ECLDATA.Qwi;    % This is cumulative  FWIT
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
    Ccomp=5e6;          % Cost to drill main bore to the top of reservoir ($)
    Cdrill=1000/0.3048;         % Drilling Cost ($/m)
    Po=80*6.2898;                % Oil price ($/SM3)
    Pg=240;      % gas price ($/MSM3)
    Pwi=8*6.2898;               % Water injected cost ($/SM3)
    Pwp=10*6.2898;               % Water Produced cost ($/SM3)
    r=0.1;                % Discount Rate
    Lwmain=2717.9016;            % Length of the main bore (m)

    % Capital Expenditure (total cost of Drill& Complete well) ($)
    Ccapex=Nwell*(Ccomp+(Lwmain).*Cdrill);

    Qgp=Qgp/(10^3); %Convert to MSM3
    Rt=Po.*Qop+Pg.*Qgp;    % Renevue  at time t
    Et=Pwp.*Qwp+Pwi.*Qwi;
    CFt=Rt-Et;               % Cash Flow at time t

    %%This is Net Present Value: (Revenue - Cost)/discount - capex
    NPV=(sum(CFt./(1+r).^t)-Ccapex);

return


