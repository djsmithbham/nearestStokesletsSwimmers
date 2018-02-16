%GENERATESPERMCONVERGENCETABLES Generates the NTABLEth table from the
%supplemental material
%
function GenerateSpermConvergenceTables(nTable)

%% Table 1 nh = 96, Qh = 96 
if nTable == 1
    
    nBeats = 1;
    f1 = 4;
    f3 = 4;
    
    f2 = [40,80,160,320,640];
    f4 = [100,200,400,800,1600];
    
    SpermConvergenceTable1 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                SpermConvergenceTable1{ii,jj} = ConvSperm(f2(ii),f1,...
                    f4(jj),f3,nBeats,0);                
            end
        end
    end
    
    save('SpermConvergenceTable1.mat','SpermConvergenceTable1','-v7.3')
end

%% Table 2 nh = 96, Qh = 600 
if nTable == 2
    
    nBeats = 1;
    f1 = 4;
    f3 = 10;
    
    f2 = [40,80,160,320,640];
    f4 = [100,200,400,800,1600];
    
    SpermConvergenceTable2 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                SpermConvergenceTable2{ii,jj} = ConvSperm(f2(ii),f1,...
                    f4(jj),f3,nBeats,0);                
            end
        end
    end
    
    save('SpermConvergenceTable2.mat','SpermConvergenceTable2','-v7.3')
end

%% Table 3 nh = 600, Qh = 600
if nTable == 3
    
    nBeats = 1;
    f1 = 10;
    f3 = 10;
    
    f2 = [40,80,160,320,640];
    f4 = [100,200,400,800,1600];
    
    SpermConvergenceTable3 = cell(length(f2),length(f4));
    
    for ii = 1:length(f2)
        for jj = 1:length(f4)
            if f4(jj) > f2(ii)
                SpermConvergenceTable3{ii,jj} = ConvSperm(f2(ii),...
                    f1,f4(jj),f3,nBeats,0);                
            end
        end
    end
    
    save('SpermConvergenceTable3.mat','SpermConvergenceTable3','-v7.3')
end

