clc
clear
close all

%  Marine Predators Algorithm (2020)
global NFE
NFE = 0;


%% Define the Problem
nVar = 10;  % number of Decision Variables
VarMin = -30; % Lower boundry
VarMax = 40;  % Upper
CostFunction = @(x)Ackley(x);  % Objective Function
%CostFunction = @(x)func_plot(F20)
%CostFunction = @(x)Sphere(x)

%% General Parameters
MaxIt = 100;  % Maximum Iteration
nPop = 40;   % Number of Population (Prey)

%% Special Parameters of MPA
FADs = 0.2;
P = 0.5;

%% initilization
fitness = zeros(nPop,1);
Prey = zeros(nPop,nVar);
for i = 1:nPop
    Prey(i,:) = unifrnd(VarMin,VarMax,[1,nVar]);
    fitness(i) = CostFunction(Prey(i,:));
end

% Find The Best
[BsetFit,iDBest] =  min(fitness);
Top_predator.Pos =  Prey(iDBest,:);
Top_predator.Fit =  BsetFit;

Xmin=repmat(ones(1,nVar).*VarMin,nPop,1);
Xmax=repmat(ones(1,nVar).*VarMax,nPop,1);

Convergence_curve=zeros(1,MaxIt); % best fitness(cost) function
stepsize=zeros(nPop,nVar);

fit_old = fitness;    Prey_old = Prey;

nfe = zeros(1,MaxIt); % number of function evaluation

%% Main Loop of MPA

for it = 1:MaxIt

    Elite = repmat(Top_predator.Pos,nPop,1);  %(Eq. 10)
    CF = (1-it/MaxIt)^(2*it/MaxIt);

    RL = 0.05*levy(nPop,nVar,1.5);   % Levy random number vector
    RB = randn(nPop,nVar);           % Brownian random number vector

    for i = 1:nPop
        for j = 1:nVar
            R = rand();
            %------------------ Phase 1 (Eq.12) -------------------
            if it<MaxIt/3
                stepsize(i,j) = RB(i,j)*(Elite(i,j)-RB(i,j)*Prey(i,j));
                Prey(i,j) = Prey(i,j) + P*R*stepsize(i,j);

                %--------------- Phase 2 (Eqs. 13 & 14)----------------
            elseif it>MaxIt/3 && it<2*MaxIt/3

                if i>nPop/2
                    stepsize(i,j) = RB(i,j)*(RB(i,j)*Elite(i,j) - Prey(i,j));
                    Prey(i,j) = Elite(i,j)+P*CF*stepsize(i,j);
                else
                    stepsize(i,j) = RL(i,j)*(Elite(i,j)-RL(i,j)* Prey(i,j));
                    Prey(i,j) = Prey(i,j)+P*R*stepsize(i,j);
                end

                %----------------- Phase 3 (Eq. 15)-------------------
            else

                stepsize(i,j) = RL(i,j)*(RL(i,j)*Elite(i,j) - Prey(i,j));
                Prey(i,j) = Elite(i,j) + P*CF*stepsize(i,j);

            end
        end
    end

    %------------------ Detecting top predator ------------------
    for i = 1:nPop

        Prey(i,:) = max(min(Prey(i,:),VarMax),VarMin);
        fitness(i) = CostFunction(Prey(i,:));

        if fitness(i)<Top_predator.Fit
            Top_predator.Fit = fitness(i);
            Top_predator.Pos = Prey(i,:);
        end
    end

    %---------------------- Marine Memory saving ----------------
    Inx = (fit_old < fitness);
    Indx = repmat(Inx,1,nVar);
    Prey = Indx.*Prey_old + ~Indx.*Prey;
    fitness = Inx.*fit_old + ~Inx.*fitness;

    fit_old = fitness;    Prey_old = Prey;

    %---------- Eddy formation and FAD's effect (Eq 16) -----------

    if rand() < FADs
        U = rand(nPop,nVar)<FADs;
        Prey = Prey + CF*((Xmin+rand(nPop,nVar).*(Xmax-Xmin)).*U);

    else
        r = rand();  Rs = size(Prey,1);
        stepsize = (FADs*(1-r)+r)*(Prey(randperm(Rs),:)-Prey(randperm(Rs),:));
        Prey = Prey + stepsize;
    end

    %------------------- Detecting top predator -----------------
    for i = 1:nPop
        Prey(i,:) = max(min(Prey(i,:),VarMax),VarMin);
        fitness(i) = CostFunction(Prey(i,:));

        if fitness(i) < Top_predator.Fit
            Top_predator.Fit = fitness(i);
            Top_predator.Pos = Prey(i,:);
        end
    end

    Convergence_curve(it) = Top_predator.Fit;

    %------------------- Marine Memory saving -------------------
    Inx = (fit_old<fitness);
    Indx = repmat(Inx,1,nVar);
    Prey = Indx.*Prey_old+~Indx.*Prey;
    fitness = Inx.*fit_old+~Inx.*fitness;

    fit_old = fitness;    Prey_old = Prey;

    nfe(it) = NFE;
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost = ' num2str(Top_predator.Fit)]);
end

%% Results
figure;
plot(1:MaxIt,Convergence_curve,'LineWidth',2);
% semilogy(nfe,BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
title('MPA Progress');


