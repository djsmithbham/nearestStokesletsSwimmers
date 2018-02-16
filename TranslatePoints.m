function X=TranslatePoints(X,a)

% performs translation on grid vector X by 3-vector a

  [X1, X2, X3]=ExtractComponents(X);
  X=[X1+a(1);X2+a(2);X3+a(3)];

