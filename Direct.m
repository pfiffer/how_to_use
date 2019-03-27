profile off
profile on -timer real
clear;
clc;
close all;
tstart = tic;

% %
% Thi code is written to find the maximum NPV just by MATLAB and Reservoir
% Simulator                                                               
% %

global NFE;
NFE = 0;

CostFunction = @(x) NPV(x);              %Cost function definition

nVar = input('Enter number of decision variables: \n');     %number of decision varables

VarSize = [1 nVar];                      % Size of Decision Variables Matrix

globalBest.Cost = -Inf;                  %Intial value for global best cost
NPV = zeros(25,25);
for i = 1:25
    for j = 1:25
        position = [i j];
        NPV(i,j) = CostFunction(position);
        if NPV(i,j) > globalBest.Cost
            globalBest.Cost = NPV(i,j);
        end  
    end
    disp(i)
end
cost = NPV;
% NPV(14:25,1:13) = flipud(NPV(1:12,1:13));
% NPV(:,14:25) = rot90(NPV(:,1:12), 2);
imagesc(rot90(flipud(NPV), -1));
xlabel('x')
ylabel('y')
cost(cost<max(cost(:))) = 0;
[m n] = find(cost);
Location = [m n];
% If model is heterogeneous
% In case of homogenous cases, comment the below lines
Perm = load('Perm-het.txt');
% Perm = ones(25,25)*20;
figure
imagesc(Perm);
hold on
plot(Location(end,1),Location(end,2),'o','MarkerSize',8,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
tElapsed = toc(tstart);

