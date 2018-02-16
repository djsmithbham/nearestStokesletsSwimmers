function x=GenerateSphereDiscr(m,a)

  % generates points on surface of sphere using cube projection
  % Input:   m,    number of subdivisions per face
  %          a,    sphere radius
  % Output:  x,    array of points

  [~,xc,yc,zc]=UniformSphereMeshGen(1,m);
  M=6*m*m;
  x=a*[reshape(permute(xc,[2,3,1]),M,1); ...
      reshape(permute(yc,[2,3,1]),M,1);reshape(permute(zc,[2,3,1]),M,1)]; 
