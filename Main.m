profile off;
profile on -timer real
clear;
clc;
close all;
tstart = tic;

%%
% Thi code is written to find the maximum NPV by PSO with limited
% search domain and Reservoir
% Simulator

global NFE;
NFE = 0;
nx = 25; ny = 25;

Perm = load('Perm-het.txt');
perm = zeros(1, nx*ny);
for i = 1:size(Perm,1)
    for j = 1:size(Perm,2)
        perm(((i-1)*25)+j) = Perm(i,j);
    end
end
Perm = zeros(ny,nx);
for  i = 1:ny
    for j =1:nx
        Perm(i,j) = perm(nx*(i-1)+j);
    end
end

%% Uncomment the below section if you want the model to be homogenous

% Perm=ones(25,25)*20;
% imagesc(Perm);
% axis image

CostFunction = @(x) NPV(x);         % Cost function definition

nVar = input('Enter number of decision variables: \n');   %number of decision varables

VarSize = [1 nVar];

VarMin = 0;     % Lower boundary of variables
VarMax = 1;     % Upper boundary of variables

%% Fast Marching Method TOF
TimeFunction = @(x, P) FMM_2(x, P);  % Time of Flight Function

position = [1 1 1 25 25 1 25 25];   % prod well position
[TOF, Index] = TimeFunction(position, Perm);
figure
imagesc(TOF);
axis image

figure
imagesc(Index)
axis image
hold on

%% We add some grids around he drainage boundary to the search domain.
%% Not only the exact drainage boundary
for i = 1:ny
    for j = 1:nx-1
        if((Index(i,j) ~= Index(i, j+1)) == 1 && (Index(i,j) ~= 0))
            Index(i,j) = 0;
            Index(i,j+1) = 0;
            Index(i,j+2) = 0;
            if(j~=1)
                Index(i,j-1) = 0;
            end
        end
    end
end

for i = 1:ny-1
    for j = 1:nx
        if((Index(i,j) ~= Index(i+1,j)) && (Index(i,j) ~= 0) && (Index(i+1,j) ~= 0 ))
            Index(i,j) = 0;
            Index(i+1,j) = 0;
            Index(i+2,j) = 0;
            if(i~=1)
                Index(i-1,j) = 0;
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
xlabel('X')
ylabel('Y')
title('Drainage boundary of four producer wells for the synthetic homogeneous case using FMM')
legend('Producer 1 drainage area' ,'Producer 2 drinage area', 'Producer 3 drainage area', 'Producer 4 drainage area')
hold on
plot(1,25, 's', 'MarkerSize', 9, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(25,1, 's', 'MarkerSize', 9, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(1,1, 's', 'MarkerSize', 9, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(25,25, 's', 'MarkerSize', 9, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
Index(Index == 0) = -1;
Index(Index ~= -1) = 0;
[x y] = find(Index);
PSO_SearchDomain = [y x]; % MATLAB indexing is different from ECLIPSE
DomainLength = length(PSO_SearchDomain);
%% PSO Parameters

MaxIt = 20;    % Maximum number of iterations

nPop = 5;      % Number of perticles

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

globalBest.Cost = -Inf;

M = 25;

% Assigning initial values for each particle
for i = 1:nPop
    
    % Initial Position
    particle(i).Position = rand(VarSize);

    % Initial Cost
    particle(i).Cost = NPV(particle(i).Position, PSO_SearchDomain);
        
    % Initial Velocity
    particle(i).Velocity = zeros(VarSize);
    
    % Updating Personal Best
    particle(i).Best.Cost = particle(i).Cost;
    particle(i).Best.Position = particle(i).Position;

    % Updating global best
    if(particle(i).Best.Cost > globalBest.Cost)
        globalBest = particle(i).Best;
    end
end

BestCost = zeros(MaxIt, 1);
Location = zeros(MaxIt, 2);
BHP = zeros(MaxIt,1);

nfe = zeros(MaxIt, 1);
tElapsed = zeros(MaxIt, 1);

%% Main loop
for It = 1:MaxIt
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
        particle(i).Cost = NPV(particle(i).Position, PSO_SearchDomain);
        
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
    BestCost(It) = globalBest.Cost;
    Location(It, :) = PSO_SearchDomain(min(floor(1+(DomainLength).*globalBest.Position(1)), DomainLength), :);
    nfe(It) = NFE;
    w = w*wdamp;
    tElapsed(It) = toc(tstart);
    
end

% now we have best position, so we can obtain the global best position's
% x and y coordinates
M = DomainLength;
globalBest.Position = PSO_SearchDomain(min(floor(1+(M).*globalBest.Position(1)), M), :);


figure
plot(nfe, BestCost, 'g--*', 'LineWidth', 2);
xlabel('NFE');
ylabel('Net Present Value - MM$');

figure
imagesc(Perm)
hold on
plot(Location(end,1),Location(end,2),'o','MarkerSize',8,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
profile viewer

