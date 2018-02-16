%GENERATEFIGURESANDTABLESFROMREPORT Generates all figures and tables from
%the manuscript
%
%% Figure 1 - Biflagellate
GenerateSwimmingFigureChlamy

%% Figure 2 - Sperm
GenerateSwimmingFigureSperm

%% Tables S1 - S5 (T2) - Chlamy convergence
for ii = 1 : 5
    GenerateChlamyConvergenceTables(ii)
end

%% Tables S6 - S8 - Sperm convergence
for ii = 1 : 3
    GenerateSpermConvergenceTables(ii)
end

%% Tables S9 - S11 - Sperm convergence with boundary
for ii = 1 : 3
    GenerateSpermBoundaryConvergenceTables(ii)
end

%% Tables S12 - S15 - Boundary point convergence with sperm
for ii = 1 : 4
    GenerateBoundaryConvergenceTables(ii)
end