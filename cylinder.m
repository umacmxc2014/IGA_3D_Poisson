addpath('./NURBS/')

knotU=[0 0 0 1 1 1];
knotV=knotU;
knotW=[0 0 1 1];
pu=2;pv=2;  pw=1;
nu=3; nv=3; nw=2;
a=sqrt(2)/2;
ndim=3;
ConPts=zeros(nu,nv,nw,ndim);
weights=zeros(nu,nv,nw);

X=[-a,0,a; -2*a   0   2*a;-a, 0, a];
Y=[ a,2*a,  a; 0,0, 0;-a, -2*a, -a];

w=[1,a,1;a,1,a;1,a,1];


for i=1:nu
   for j=1:nv
ConPts(i,j,1,:)=[X(i,j),Y(i,j),0];
weights(i,j,1)=w(i,j);
end
end

for i=1:nu
   for j=1:nv
ConPts(i,j,2,:)=[X(i,j),Y(i,j),1];
weights(i,j,2)=w(i,j);
end
end

% S=PointOnBspVolume(ConPts,U,pu,0.1,V,pv,0.2,W,pw,0.3);

 [S,DF,W,DW]=NurbsVolume(ConPts,weights,knotU,pu,0.1,knotV,pv,0.,knotW,pw,0.3);

 t=2;
 
 [Q,wbar,knotUbar,knotVbar,knotWbar]=IGADegreeElevVolume(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,t);


 [S1,DF1,W1,DW1]=NurbsVolume(Q,wbar,knotUbar,pu+t,0.1,knotVbar,pv+t,0.2,knotWbar,pw+t,0.3);


 [W,DW, D2W,F,DF,D2F]=NurbsVolumeDers(ConPts,knotU,knotV,knotW,weights,pu,0.1,pv,0.2,pw,0.3);