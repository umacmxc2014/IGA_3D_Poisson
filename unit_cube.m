knotU=[0 0   1 1];
knotV=knotU;
knotW=[0 0 1 1];
pu=1;pv=1;pw=1;
nu=2; nv=2; nw=2;

ndim=3;
ConPts=zeros(nu,nv,nw,ndim);
weights=zeros(nu,nv,nw);

X=[0 0;1 1];
Y=[0 1;0 1];

w=[1 1;1 1];


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

 [S,DF,W,DW]=NurbsVolume(ConPts,weights,knotU,pu,0.1,knotV,pv,0.2,knotW,pw,0.3)