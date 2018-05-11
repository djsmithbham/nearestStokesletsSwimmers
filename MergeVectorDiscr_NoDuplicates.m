function X=MergeVectorDiscr_NoDuplicates(Xa,Xb)

% merges two discretisations made from column vectors Xa, Xb, assuming of the form [xcoords;ycoords;zcoords]
% uses 'union' to prevent duplication

[Xa1,Xa2,Xa3]=ExtractComponents(Xa);
[Xb1,Xb2,Xb3]=ExtractComponents(Xb);

Xtemp=union([Xa1 Xa2 Xa3],[Xb1 Xb2 Xb3],'rows');
X=[Xtemp(:,1);Xtemp(:,2);Xtemp(:,3)];


