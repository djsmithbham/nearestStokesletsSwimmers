function [s,xm,ym,zm]=UniformSphereMeshGen(a,n)

s=linspace(-1+1/n,1-1/n,n)';

c=ones(n,n);
xi=s*ones(1,n);
et=ones(n,1)*s';

xm=zeros(6,n,n);
ym=zeros(6,n,n);
zm=zeros(6,n,n);

xm(1,:,:)=c;
ym(1,:,:)=xi;
zm(1,:,:)=et;

xm(2,:,:)=-xi;
ym(2,:,:)=c;
zm(2,:,:)=et;

xm(3,:,:)=xi;
ym(3,:,:)=et;
zm(3,:,:)=c;

xm(4,:,:)=xi;
ym(4,:,:)=-et;
zm(4,:,:)=-c;

xm(5,:,:)=xi;
ym(5,:,:)=-c;
zm(5,:,:)=et;

xm(6,:,:)=-c;
ym(6,:,:)=-xi;
zm(6,:,:)=et;

[xm, ym, zm]=SphereProjector(xm,ym,zm,a);

%figure(1);hold on;plot3(reshape(xm,[],1),reshape(ym,[],1),reshape(zm,[],1),'.');
