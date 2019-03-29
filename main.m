clear;
clc
tic

% If random initial well place is needed
x = [0 0 0 1 1 0 1 1];

[DelpT, DelpReal, TOF] = FMM(x);
%%
% TOF is the diffusive Time of Flight %%
%%

figure
imagesc(TOF)
toc