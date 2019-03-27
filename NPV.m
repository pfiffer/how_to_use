function NPV  = NPV(x)
    global NFE
    % NFE=0;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1
	
	k=x;

	ECLDATA.CompLoc(1:2)=k;
	ECLDATA.CompLoc(3:4)=1;
	ECLDATA.WHLoc=k;
	ECLDATA=Eclipse(ECLDATA);

	Qopcum=ECLDATA.Qop;    % This is FGPT
	Qwpcum=ECLDATA.Qwp;    % This FWPT
	Qwicum=ECLDATA.Qwi;    % This is cumulative  FWIT
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
	Nwell=numel(x)/2; % number of wells (4production wells)
	% Economic Parameters :
	Cdrill=1000;          % Drilling Cost ($/ft)
	Ccomp = 5*10^6;		  % Completeion costs ($)
	Po=80;                % Oil price ($/STB)
	Pg=240*0.3048^3;      % gas price ($/MScf)
	Pwi=8;                % Water injected cost ($/STB)
	Pwp=10;               % Water Produced cost ($/STB)
	r=0.1;                % Discount Rate
	Lwmain=8917;          % Length of the main bore 0.3048=ft 

	% Capital Expenditure (total cost of Drill& Complete well) ($)
	Ccapex=Nwell*(Ccomp + (Lwmain)*Cdrill);

	Qgp=0;
	Rt=Po.*Qop+Pg.*Qgp;      % Renevue  at time t
	Et=Pwp.*Qwp+Pwi.*Qwi;
	CFt=Rt-Et;               % Cash Flow at time t

	%%This is Net Present Value: (Revenue - Cost)/discount - capex
	NPV=sum(CFt./(1+r).^t)-Ccapex;

return
