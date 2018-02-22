%GENERATECHLAMYCONVERGENCETABLES Generates the NTABLEth table from 
%the supplemental material
%
function GenerateChlamyConvergenceTables(nTable)

%% Table 1 nh = 96, Qh = 600
if nTable == 1
    
    nBeats = 1;
    f1 = 4;
    f3 = 10;
    
    f2 = [40,80,160,320];
    f4 = [100,200,400,800,1600];
    
    ChlamyConvergenceTable1 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                ChlamyConvergenceTable1{ii,jj} = ConvChlamy(f2(ii),f1, ...
                    f4(jj),f3,nBeats);                
            end
        end
    end
    
    save('ChlamyConvergenceTable1.mat','ChlamyConvergenceTable1','-v7.3')
end

%% Table 2 nh = 600, Qh = 600
if nTable == 2
    
    nBeats = 1;
    f1 = 10;
    f3 = 10;
    
    f2 = [40,80,160,320];
    f4 = [100,200,400,800,1600];
    
    ChlamyConvergenceTable2 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                ChlamyConvergenceTable2{ii,jj} = ConvChlamy(f2(ii),f1, ...
                    f4(jj),f3,nBeats);                
            end
        end
    end
    
    save('ChlamyConvergenceTable2.mat','ChlamyConvergenceTable2','-v7.3')
end

%% Table 3 nh = 96, Qh = 2646
if nTable == 3
    
    nBeats = 1;
    f1 = 4;
    f3 = 21;
    
    f2 = [40,80,160,320];
    f4 = [100,200,400,800,1600];
    
    ChlamyConvergenceTable3 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                ChlamyConvergenceTable3{ii,jj} = ConvChlamy(f2(ii),f1, ...
                    f4(jj),f3,nBeats);              
            end
        end
    end
    
    save('ChlamyConvergenceTable3.mat','ChlamyConvergenceTable3','-v7.3')
end

%% Table 4 Nystrom
if nTable == 4
    
    nBeats = 1;
    f1 = [4,10,21];
    f2 = [40,100,200,400];
    
    ChlamyConvergenceTable4 = cell(length(f1),length(f2));
    
    for ii = 1:length(f1)
        for jj = 1:length(f2)
                ChlamyConvergenceTable4{ii,jj} = ConvChlamy(f2(ii),f1(jj), ...
                    f1(ii),f2(jj),nBeats);              
        end
    end
    
    save('ChlamyConvergenceTable4.mat','ChlamyConvergenceTable4','-v7.3')
end

