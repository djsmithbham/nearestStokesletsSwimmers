%GENERATEBOUNDARYCONVERGENCETABLES Generates the NTABLEth table from 
%the supplemental material
%
function GenerateBoundaryConvergenceTables(nTable)

%% Table 1 nx = 16, ny = 15
if nTable == 1
    
    nBeats = 1;
    f1 = 16;
    f2 = 15;
    
    f3 = [16,32,64,128,256];
    f4 = [15,30,60,120,240];
    
    BoundaryConvergenceTable1 = cell(length(f3),length(f4));
    
    for ii = 1:length(f3)
        for jj = 1:length(f4)
                BoundaryConvergenceTable1{ii,jj} ...
                    = ConvBoundary(f1,f2,f3(ii),f4(jj),nBeats);                
        end
    end
    
    save('BoundaryConvergenceTable1.mat','BoundaryConvergenceTable1', ...
        '-v7.3')
end
%% Table 2 nx = 32, ny = 15
if nTable == 2
    
    nBeats = 1;
    f1 = 32;
    f2 = 15;
    
    f3 = [32,64,128,256];
    f4 = [15,30,60,120,240];
    
    BoundaryConvergenceTable2 = cell(length(f3),length(f4));
    
    for ii = 1:length(f3)
        for jj = 1:length(f4)
                BoundaryConvergenceTable2{ii,jj} ...
                    = ConvBoundary(f1,f2,f3(ii),f4(jj),nBeats);                  
        end
    end
    
    save('BoundaryConvergenceTable2.mat','BoundaryConvergenceTable2', ...
        '-v7.3')
end
%% Table 3 nx = 16, ny = 30
if nTable == 3
    
    nBeats = 1;
    f1 = 16;
    f2 = 30;
    
    f3 = [16,32,64,128,256];
    f4 = [30,60,120,240];
    
    BoundaryConvergenceTable3 = cell(length(f3),length(f4));
    
    for ii = 1:length(f3)
        for jj = 1:length(f4)
                BoundaryConvergenceTable3{ii,jj} ...
                    = ConvBoundary(f1,f2,f3(ii),f4(jj),nBeats);                                
        end
    end
    
    save('BoundaryConvergenceTable3.mat','BoundaryConvergenceTable3', ...
        '-v7.3')
end
%% Table 4 nx = 32, ny = 30
if nTable == 4
    
    nBeats = 1;
    f1 = 32;
    f2 = 30;
    
    f3 = [32,64,128,256];
    f4 = [30,60,120,240];
    
    BoundaryConvergenceTable4 = cell(length(f3),length(f4));
    
    for ii = 1:length(f3)
        for jj = 1:length(f4)
                BoundaryConvergenceTable4{ii,jj} ...
                    = ConvBoundary(f1,f2,f3(ii),f4(jj),nBeats);                                
        end
    end
    
    save('BoundaryConvergenceTable4.mat','BoundaryConvergenceTable4', ...
        '-v7.3')
end