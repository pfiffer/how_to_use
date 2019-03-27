profile off
profile on -timer real
clear;
clc;
close all;
tstart = tic;

%%
% Thi code is written to find the maximum NPV just by PSO and Reservoir
% Simulator                                                              

global NFE;
NFE = 0;
load('initRand.mat');
%CostFunction = @(x) NPV(x);              % Cost function definition

nVar = input('Enter number of decision variables: \n');     % number of decision varables

VarSize = [1 nVar];                      % Size of Decision Variables Matrix
nx = 19; ny=28;
Perm = load('Permx.txt');
Perm = reshape(Perm, [nx ny 5]); 
% perm = zeros(1, nx*ny);
% for i = 1:size(Perm,1)
%     for j = 1:size(Perm,2)
%         perm(((i-1)*1)+j) = Perm(i,j);
%     end
% end
% Perm = zeros(ny,nx);
% for  i = 1:ny
%     for j =1:nx
%         Perm(i,j) = perm(nx*(i-1)+j);
%     end
% end

nwp = input('Enter number of production wells: ');  %number of production wells

VarMax = 1;                              % Upper boundary of variables
VarMin = 0;                              % Lower boundary of variables

%% PSO Parameters
MaxIt = 20;         % Maximum number of iterations

nPop = 10;          % Number of particles


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

%% Initialization of parameters
empty_particle.Position = [];
empty_particle.Cost = [];
empty_particle.Velocity = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

particle = repmat(empty_particle, nPop, 1);

% These for every PSO Iteraion
BestCost = zeros(MaxIt+1, 1);
Location = zeros(MaxIt+1, nVar);
BHP = zeros(MaxIt+1, 1);
NetPresentValue = zeros(MaxIt+1, 1);

nfe = zeros(MaxIt+1, 1);
tElapsed = zeros(MaxIt+1, 1);

globalBest.Cost = -Inf;                  % Intial value for global best cost
%% Assigning initial values
for i = 1:nPop
    %initial Position
    particle(i).Position = rand(VarSize); %Randd(i,1:4);
    
    %Evaluation
    particle(i).Cost = NPV(particle(i).Position, nx, ny, nwp);
    
    %initial velocity
    particle(i).Velocity = zeros(VarSize);
    
    %Updating personal best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    %Updating global best
    if(particle(i).Best.Cost > globalBest.Cost)
        globalBest = particle(i).Best;
    end
    
    BestCost(1) = globalBest.Cost;
    if(nVar > 2)
        for j=1:nVar/2
            if(mod(j,2) ~= 0)
                Location(1, j:2:end) = min(floor(1+(nx).*globalBest.Position(j:2:end)), nx);
            else
                Location(1, j:2:end) = min(floor(1+(ny).*globalBest.Position(j:2:end)), ny);
            end
        end
    else
        Location(1, 1:2:end) = min(floor(1+(nx).*globalBest.Position(1:2:end)), nx);
        Location(1, 1:2:end) = min(floor(1+(nx).*globalBest.Position(1:2:end)), nx);
    end
    nfe(1) = NFE;
    tElapsed(1) = toc(tstart);
    NetPresentValue(1) = BestCost(1);
end



%% PSO main lopp

for it = 1:MaxIt
    for i = 1:nPop
        % Update velocity
        particle(i).Velocity = w.*particle(i).Velocity...
        + c1.*rand(VarSize).*(particle(i).Best.Position-particle(i).Position)...
        + c2.*rand(VarSize).*(globalBest.Position-particle(i).Position);
        
        % Applying velocity limits
        particle(i).Velocity = max(particle(i).Velocity, VelMin);
        particle(i).Velocity = min(particle(i).Velocity, VelMax);
    
        %update position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        %velocity mirror effect
        IsOutside = (particle(i).Position < VarMin | particle(i).Position > VarMax);
        particle(i).Velocity(IsOutside) = -particle(i).Velocity(IsOutside);
        
        % Applying position limits
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);
        
        % Evaluating cost with new position
        particle(i).Cost = NPV(particle(i).Position, nx, ny, nwp);
        
        % Update personal best
        if(particle(i).Cost >= particle(i).Best.Cost)
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.Position = particle(i).Position;
            
            % Update Global Best
            if(particle(i).Best.Cost > globalBest.Cost)
                globalBest = particle(i).Best;
            end
        end
    end
    BestCost(it+1) = globalBest.Cost;
    if(nVar > 2)
        for i=1:nVar/2
            if(mod(i,2) ~= 0)
                Location(it+1, i:2:end) = min(floor(1+(nx).*globalBest.Position(i:2:end)), nx);
            else
                Location(it+1, i:2:end) = min(floor(1+(ny).*globalBest.Position(i:2:end)), ny);
            end
        end
    else
        Location(it+1, 1:2:end) = min(floor(1+(nx).*globalBest.Position(1:2:end)), nx);
        Location(it+1, 1:2:end) = min(floor(1+(nx).*globalBest.Position(1:2:end)), nx);
    end
    nfe(it+1) = NFE;
    tElapsed(it+1) = toc(tstart);
    NetPresentValue(it+1) = BestCost(it+1);
    
    %BHP(it)=FMMNew( GlobalBest.Position);
    w = w*wdamp;
    tElapsed(it+1) = toc(tstart);
end

% now we have best position, so we can obtain the global best position's
% x and y coordinates
globalBest.Position = Location(end, :);

figure
plot(nfe, BestCost, 'g--*', 'LineWidth', 2);
xlabel('NFE');
ylabel('Net Present Value - $');

figure
imagesc((Perm(:,:,1)'));
axis image
hold on
for i=1:nVar/2
    plot(Location(end,2*i-1),Location(end,2*i),'o','MarkerSize',8,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
end

profile viewer
