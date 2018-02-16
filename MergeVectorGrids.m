function X=MergeVectorGrids(Xa,Xb)

% merges two grids made from column vectors X1, X2, 
%      assuming of the form [xcoords;ycoords;zcoords]

[Xa1, Xa2, Xa3]=ExtractComponents(Xa);
[Xb1, Xb2, Xb3]=ExtractComponents(Xb);

X=[Xa1;Xb1;Xa2;Xb2;Xa3;Xb3];
