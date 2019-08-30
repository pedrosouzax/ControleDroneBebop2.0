%clc, close all, clear all

FitFcn = @Main;
nvars = 3;


options = optimoptions('ga','FunctionTolerance', 4, 'MaxGenerations', 5, 'MaxStallGenerations', 5, 'PopulationSize', 10 )

% opt = gaoptimset('Generations',3);
 
[x,fval,exitflag, output] = ga(FitFcn, nvars,[],[],[],[],[],[],[],options)