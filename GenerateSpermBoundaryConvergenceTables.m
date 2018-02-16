%GENERATESPERMBOUNDARYCONVERGENCETABLES Generates the NTABLEth table from 
%the supplemental material
%
function GenerateSpermBoundaryConvergenceTables(nTable)

%% Table 1 nh = 96, Qh = 96 
if nTable == 1
    
    nBeats = 1;
    f1 = 4;
    f3 = 4;
    
    f2 = [40,80,160,320,640];
    f4 = [100,200,400,800,1600];
    
    SpermBoundaryConvergenceTable1 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                SpermBoundaryConvergenceTable1{ii,jj} = ConvSperm(f2(ii),...
                    f1,f4(jj),f3,nBeats,1);                
            end
        end
    end
    
    save('SpermBoundaryConvergenceTable1.mat', ...
        'SpermBoundaryConvergenceTable1','-v7.3')
end

%% Table 2 nh = 96, Qh = 600 
if nTable == 2
    
    nBeats = 1;
    f1 = 4;
    f3 = 10;
    
    f2 = [40,80,160,320,640];
    f4 = [100,200,400,800,1600];
    
    SpermBoundaryConvergenceTable2 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                SpermBoundaryConvergenceTable2{ii,jj} = ConvSperm( ...
                    f2(ii),f1,f4(jj),f3,nBeats,1);                
            end
        end
    end
    
    save('SpermBoundaryConvergenceTable2.mat', ...
        'SpermBoundaryConvergenceTable2','-v7.3')
end

%% Table 3 nh = 600, Qh = 600
if nTable == 3
    
    nBeats = 1;
    f1 = 10;
    f3 = 10;
    
    f2 = [40,80,160,320,640];
    f4 = [100,200,400,800,1600];
    
    SpermBoundaryConvergenceTable3 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                SpermBoundaryConvergenceTable3{ii,jj} = ConvSperm( ...
                    f2(ii),f1,f4(jj),f3,nBeats,1);                
            end
        end
    end
    
    save('SpermBoundaryConvergenceTable3.mat', ...
        'SpermBoundaryConvergenceTable3','-v7.3')
end

