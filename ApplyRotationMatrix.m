function x=ApplyRotationMatrix(M,xi)

[xi1,xi2,xi3]=ExtractComponents(xi);
y=M*[xi1';xi2';xi3'];
x=[y(1,:)';y(2,:)';y(3,:)'];
