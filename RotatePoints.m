function X=RotatePoints(X,X0,theta,axx)
  [X1 X2 X3]=ExtractComponents(X);
  [O1 O2 O3]=ExtractComponents(X0);
  R=RotationMatrix(theta,axx);
  Ze=0*X1';
  X=R*[(X1-O1)';(X2-O2)';(X3-O3)']+[O1+Ze;O2+Ze;O3+Ze];
  X=reshape(X',[],1);

