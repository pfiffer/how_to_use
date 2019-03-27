profile off;
profile on -timer real
clear;
clc;
close all;
tstart = tic;

%%
% Thi code is written to find the maximum NPV using PSO-FMM with limited
% search domain and Reservoir Simulator

global NFE;
NFE = 0;

nx = 60; ny = 220;
Perm = load('Perm-het.txt');
perm = zeros(1, nx*ny);
for i = 1:size(Perm,1)
    for j = 1:size(Perm,2)
        perm(((i-1)*6)+j) = Perm(i,j);
    end
end
Perm = zeros(ny,nx);
for  i = 1:ny
    for j =1:nx
        Perm(i,j) = perm(nx*(i-1)+j);
    end
end

 CostFunction = @(x) NPV(x);         % Cost function definition

 nVar = input('Enter number of decision variables: \n');   %number of decision varables

 VarSize = [1 nVar];

 VarMin = 0;     % Lower boundary of variables
 VarMax = 1;     % Upper boundary of variables

%% Fast Marching Method TOF
 TimeFunction = @(x, Perm) FMM_2(x, Perm);  % Time of Flight Function

position = [1 1 1 floor(ny/2) 1 ny nx 1 nx floor(ny/2) nx ny];   % prod well position
nwp = size(position,2)/2;       % number of prod wells
[TOF, Index] = TimeFunction(position, Perm);
TOF = real(TOF);
figure
imagesc((TOF));
axis image
figure
imagesc(Index)
axis image

%% We add some grids around he drainage boundary to the search domain.
%% Not only the exact drainage boundary
for i = 1:ny
    for j = 1:nx-4
        if((Index(i,j) ~= Index(i, j+1)) == 1 && (Index(i,j) ~= 0) && Index(i,j+2) ~= 0)
            Index(i,j) = 0;
            Index(i,j+1:j+4) = 0;
            if(j~=1)
                Index(i,j-1) = 0;
            end
            % More the number of total grids, more the number neighbors
            if(j-2>2)
                Index(i,j-2) = 0;
                Index(i,j-3) = 0;
            end
        end
    end
end

for i = 1:ny-3
    for j = 1:nx
        if((Index(i,j) ~= Index(i+1,j)) && (Index(i,j) ~= 0) && (Index(i+2,j) ~= 0 ))
            Index(i,j) = 0;
            Index(i+1:i+4,j) = 0;
            if(i~=1)
                Index(i-1,j) = 0;
            end
            % More the number of total grids, more the number neighbors
            if(i-2>2)
                Index(i-2,j) = 0;
                Index(i-3,j) = 0;
            end
        end
    end
end
for i = 1:ny-1
    for j = 2:nx-1
        if(Index(i,j-1) == 0 && Index(i,j+1) == 0 && Index(i,j) ~= 0) %% || ((Index(i-1,j) == Index(i+1,j)) && i > 0)))
            Index(i,j) = 0;
        end
    end
end

figure
imagesc(Index)
axis image
Index(Index == 0) = -1;
Index(Index ~= -1) = 0;
[x y] = find(Index);
PSO_SearchDomain = [y x];       % Because MATLBAD indexing is different form ECLIPSE
DomainLength = length(PSO_SearchDomain);
%% PSO Parameters

MaxIt = 20;     % Maximum number of iterations
nPop = 15;      % Number of perticles

% Constriction Coefficients
phi1 = 2.05;
phi2 = 2.05;
phi = phi1+phi2;
chi = 2/(phi-2+sqrt(phi^2-4*phi));
w = chi;          % Inertia Weight
wdamp = 1;        % Inertia Weight Damping Ratio
c1 = chi*phi1;    % Personal Learning Coefficient
c2 = chi*phi2;    % Global Learning Coefficient

% Velocity limits
VelMax = 0.1*(VarMax-VarMin);
VelMin = -VelMax;

% Initialization of parameters
empty_partcle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];
particle = repmat(empty_particle, nPop, 1);

% These for every PSO Iteraion
BestCost = zeros(MaxIt+1, 1);
Location = zeros(MaxIt+1, nVar);
BHP = zeros(MaxIt+1, 1);

nfe = zeros(MaxIt+1, 1);
tElapsed = zeros(MaxIt+1, 1);
NetPresentValue = zeros(MaxIt+1, 1);

globalBest.Cost = -Inf;

% Assigning initial values for each particle
for i = 1:nPop
    
    % Initial Position
    particle(i).Position = rand(VarSize);

    % Initial Cost
    particle(i).Cost = NPV(particle(i).Position, PSO_SearchDomain, nwp);
        
    % Initial Velocity
    particle(i).Velocity = zeros(VarSize);
    
    % Updating Personal Best
    particle(i).Best.Cost = particle(i).Cost;
    particle(i).Best.Position = particle(i).Position;

    % Updating global best
    if(particle(i).Best.Cost > globalBest.Cost)
        globalBest = particle(i).Best;
    end
     % including inital values in the output
    BestCost(1) = globalBest.Cost;
    for ii = 1:nVar/2
         Location(1, 2*ii-1:2*ii) = PSO_SearchDomain(min(floor(1+(DomainLength).*globalBest.Position(2*ii-1)),(DomainLength)), :);
    end
    nfe(1) = NFE;
    tElapsed(1)=toc(tstart);
    NetPresentValue(1) = BestCost(1);
end
figure
plot(nfe(1), BestCost(1), 'g--*', 'LineWidth', 2);
xlabel('NFE');
ylabel('Net Present Value - $');
hold on

%% Main loop
for it = 20:MaxIt
    for i = 1:nPop
        
        % Updating Velocity
        particle(i).Velocity = w.*particle(i).Velocity...
        + c1.*rand(VarSize).*(particle(i).Best.Position-particle(i).Position)...
        + c2.*rand(VarSize).*(globalBest.Position-particle(i).Position);
    
        % Applying velocity limits
        particle(i).Velocity = max(particle(i).Velocity, VelMin);
        particle(i).Velocity = min(particle(i).Velocity, VelMax);
        
        % Updating Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity mirror effect
        IsOutside = (particle(i).Position > VarMax | particle(i).Position < VarMin);
        particle(i).Velocity(IsOutside) = -particle(i).Velocity(IsOutside);
        
        % Applying position limits
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);
               
        % Evaluating cost
        particle(i).Cost = NPV(particle(i).Position, PSO_SearchDomain, nwp);
        
        % Updating Personal best
        if(particle(i).Cost > particle(i).Best.Cost)
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.Position = particle(i).Position;
            
            % Update global best
            if(particle(i).Best.Cost > globalBest.Cost)
                globalBest = particle(i).Best;
            end
        end
    end
    BestCost(it+1) = globalBest.Cost;
%     M = 25;
%     N = 50;
    for ii = 1:nVar/2
         Location(it+1, 2*ii-1:2*ii) = PSO_SearchDomain(min(floor(1+(DomainLength).*globalBest.Position(2*ii-1)),(DomainLength)), :);
    end
    nfe(it+1) = NFE;
    
    %BHP(it)=FMMNew( GlobalBest.Position);
    w = w*wdamp;
    tElapsed(it+1) = toc(tstart);
    NetPresentValue(it+1) = BestCost(it+1);
    save('Result_PSOFMM_2Inj_BHP_2(Extended).mat', 'particle', 'globalBest', 'nfe', 'BestCost', 'tElapsed', 'Location');

    plot(nfe(it+1), BestCost(it+1), 'g--*', 'LineWidth', 2);
    xlabel('NFE');
    ylabel('Net Present Value - $');
    hold on
    
end

% now we have best position, so we can obtain the global best position's
% x and y coordinates
globalBest.Position = Location(end, :);

figure
plot(nfe, BestCost, 'g--*', 'LineWidth', 2);
xlabel('NFE');
ylabel('Net Present Value - $');

figure
imagesc(log10(Perm));
axis image
hold on
for i=1:nVar/2
    plot(Location(end,2*i-1),Location(end,2*i),'s','MarkerSize',9,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
end
profile viewer

