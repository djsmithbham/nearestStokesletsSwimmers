function TablesForReport_SphereAugmentedDiscr(epsilon,MColl,QColl)

%Stokes' Law test case: this version with collpts and quadpts different and force points difference

%force discretised via nearest neighbour interpolation


U=[1,0,0];
Om=[1,0,0];
X0=[0,0,0];

a=1;
domain='i';
blockSize=0.2;

for iM=1:length(MColl)
    for iQ=1:length(QColl)
        tic;
        fprintf('iM = %i, iQ = %i\n',iM,iQ);
        [x,X]=GenerateSpherePoints(MColl(iM),QColl(iQ),a);
        
        % augment quadrature points with force points
        X=MergeVectorDiscr_NoDuplicates(x,X);
        
        [FTr,MTr,fTr,condno]=SolveRigidResistance(x,X,X0,U,     [0,0,0],epsilon,domain,blockSize);
        [FRo,MRo,fRo,~]     =SolveRigidResistance(x,X,X0,[0,0,0],Om,    epsilon,domain,blockSize);
        FTrError(iM,iQ)=norm(FTr-[6*pi*a  ;0;0])/6/pi/a;
        MRoError(iM,iQ)=norm(MRo-[8*pi*a^3;0;0])/8/pi/a^3;
        hf(iM)=CalcDiscr_h(x,blockSize);
        hq(iQ)=CalcDiscr_h(X,blockSize);
        deltaClosest(iM,iQ)=FindSmallestDistanceBetweenDiscretizations(X,x,blockSize);
        Dof(iM)=length(x);
        NQ(iQ)=length(X)/3;
        NM(iM)=length(x)/3;
        walltime(iM,iQ)=toc;
        condNo(iM,iQ)=condno;
    end
end

fid=fopen(['condno_transrot_test_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'    \\begin{tabular}{ccc|');
for iQ=1:length(QColl)
    fprintf(fid,'c');
end
fprintf(fid,'}\n    ');
fprintf(fid,'          &     &    \\(Q   \\)');fprintf(fid,'& %i',NQ);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'\\(N  \\) & DOF &           \\\\[0.3em]  \\hline \n    ');
fprintf(fid,'          &     &           \\\\[-0.8em] \n    ');
for iM=1:length(MColl)
    fprintf(fid,'%i & %i &  ',NM(iM),Dof(iM));
    fprintf(fid,'& %.3f',condNo(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid)

fid=fopen(['walltime_transrot_test_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'    \\begin{tabular}{ccc|');
for iQ=1:length(QColl)
    fprintf(fid,'c');
end
fprintf(fid,'}\n    ');
fprintf(fid,'          &     &    \\(Q   \\)');fprintf(fid,'& %i',NQ);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'\\(N  \\) & DOF &           \\\\[0.3em]  \\hline \n    ');
fprintf(fid,'          &     &           \\\\[-0.8em] \n    ');
for iM=1:length(MColl)
    fprintf(fid,'%i & %i &  ',NM(iM),Dof(iM));
    fprintf(fid,'& %.3f',walltime(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid)


fid=fopen(['translation_test_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'%%Translating sphere, infinite fluid test\n');
fprintf(fid,'%%Radius, %f\n',a);
fprintf(fid,'%%Regularisation parameter, %f\n',epsilon);
fprintf(fid,'    \\begin{tabular}{cccc|');
for iQ=1:length(QColl)
    fprintf(fid,'c');
end
fprintf(fid,'}\n    ');
fprintf(fid,'          &     &           & \\(Q   \\)');fprintf(fid,'& %i',NQ);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'          &     &           & \\(h_q \\)');fprintf(fid,'& %1.4f',hq);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'\\(N  \\) & DOF & \\(h_f\\) &          \\\\[0.3em]  \\hline \n    ');
fprintf(fid,'          &     &           &          \\\\[-0.8em] \n    ');
for iM=1:length(MColl)
    fprintf(fid,'%i & %i & %1.4f & ',NM(iM),Dof(iM),hf(iM));
    fprintf(fid,'& %1.4f',FTrError(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid)

fid=fopen(['closest_point_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'%%Translating sphere, infinite fluid test\n');
fprintf(fid,'%%Radius, %f\n',a);
fprintf(fid,'%%Regularisation parameter, %f\n',epsilon);
fprintf(fid,'    \\begin{tabular}{cccc|');
for iQ=1:length(QColl)
    fprintf(fid,'c');
end
fprintf(fid,'}\n    ');
fprintf(fid,'          &     &           & \\(Q   \\)');fprintf(fid,'& %i',NQ);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'          &     &           & \\(h_q \\)');fprintf(fid,'& %1.4f',hq);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'\\(N  \\) & DOF & \\(h_f\\) &          \\\\[0.3em]  \\hline \n    ');
fprintf(fid,'          &     &           &          \\\\[-0.8em] \n    ');
for iM=1:length(MColl)
    fprintf(fid,'%i & %i & %1.4f & ',NM(iM),Dof(iM),hf(iM));
    fprintf(fid,'& %1.4f',deltaClosest(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid)

fid=fopen(['rotation_test_' num2str(epsilon) '.txt'],'w');
fprintf(fid,'%%Rotating sphere, infinite fluid test\n');
fprintf(fid,'%%Radius, %f\n',a);
fprintf(fid,'%%Regularisation parameter, %f\n',epsilon);
fprintf(fid,'    \\begin{tabular}{cccc|');
for iQ=1:length(QColl)
    fprintf(fid,'c');
end
fprintf(fid,'}\n    ');
fprintf(fid,'          &     &           & \\(Q   \\)');fprintf(fid,'& %i',NQ);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'          &     &           & \\(h_q \\)');fprintf(fid,'& %1.4f',hq);fprintf(fid,'\\\\ \n    ');
fprintf(fid,'\\(N  \\) & DOF & \\(h_f\\) &          \\\\[0.3em]  \\hline \n    ');
fprintf(fid,'          &     &           &          \\\\[-0.8em] \n    ');
for iM=1:length(MColl)
    fprintf(fid,'%i & %i & %1.4f & ',NM(iM),Dof(iM),hf(iM));
    fprintf(fid,'& %1.4f',MRoError(iM,:));
    fprintf(fid,'\\\\ \n    ');
end
fprintf(fid,'\\end{tabular}\n');
fclose(fid);


