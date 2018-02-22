%GENERATEFIGURESANDTABLESFROMREPORT Generates all figures and tables from
%the manuscript
%
%% Figure 1 - Biflagellate
GenerateSwimmingFigureChlamy

%% Figure 2 - Sperm
GenerateSwimmingFigureSperm

%% Tables S1 - S3, T2,T3 - Chlamy convergence
for ii = 1 : 4
    GenerateChlamyConvergenceTables(ii)
end

%% Tables S4 - S6 - Sperm convergence
for ii = 1 : 3
    GenerateSpermConvergenceTables(ii)
end

%% Tables S7 - S9 - Sperm convergence with boundary
for ii = 1 : 3
    GenerateSpermBoundaryConvergenceTables(ii)
end

%% Tables S10 - S13 - Boundary point convergence with sperm
for ii = 1 : 4
    GenerateBoundaryConvergenceTables(ii)
end