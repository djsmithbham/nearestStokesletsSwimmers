function R=RotationMatrix(theta,axx)

Theta=theta*ones(1,1);
I=ones(1,1);
Z=zeros(1,1);

switch axx
       case 1
       	    R=[I Z Z; Z cos(Theta) -sin(Theta); Z sin(Theta) cos(Theta)];
       case 2
       	    R=[cos(Theta) Z sin(Theta); Z I Z; -sin(Theta) Z cos(Theta)];
       case 3
       	    R=[cos(Theta) -sin(Theta) Z; sin(Theta) cos(Theta) Z; Z Z I];
end
