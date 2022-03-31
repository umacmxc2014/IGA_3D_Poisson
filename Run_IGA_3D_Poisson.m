 % test_case = 'rectangle';

  test_case = 'cylinder';

if strcmp(test_case,'rectangle')
  
%% Case I: rectangle domain
nu=2; nv=2; nw=2;
pu=1; pv=1; pw=1;
 knotU=[0 0  1 1];knotV=[0 0  1 1];knotW=[0 0 1 1];
 
 t=1;  % The (pu+t) is the ultimate degree of B splines basis functions

 DIM=3;

X = [0 0;1 1];
Y = [0 1;0 1];
ConPts=zeros(nu,nv,nw,DIM);

for i=1:nu
    for j=1:nv
        ConPts(i,j,1,1)=X(i,j);
        ConPts(i,j,2,1)=X(i,j);
        ConPts(i,j,1,2)=Y(i,j);
        ConPts(i,j,2,2)=Y(i,j);
        ConPts(i,j,1,3)=0;
        ConPts(i,j,2,3)=1;
    end
end

weights = ones(nu,nv,nw);

Refinement=[1,2,3,4];

n_refine=length(Refinement);
err=zeros(n_refine,2);
n_dofs=zeros(n_refine,1);



for i=1:size(Refinement,2)

   [err(i,:),n_dofs(i)]=IGA_3D_Poisson(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,Refinement(i),t, test_case);
 
 %  [err(i,:),n_dofs(i)]=IGA_3D_Poisson_Store_Basis(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,Refinement(i),t, test_case);

end
format short e
disp('The degree of NURBS basis is ')
disp(pu+t)
 disp('================================================================')
 disp(['The L2 norm error            ||       The H1 norm error'  ])
 
 disp('----------------------------------------------------------------')
 disp([err])
 disp('----------------------------------------------------------------')
disp('The #DOFs:')
disp(n_dofs)
disp('-----------------------------------------------------------------')
  disp('The convergence order in L2, H1 are:  ')
 
  disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2)])
 disp('================================================================')



%%

end

%% Case I: quarter domain

if strcmp(test_case,'cylinder')

 t=0;
       
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

Refinement=[1,2,3,4];

n_refine=length(Refinement);
err=zeros(n_refine,2);
n_dofs=zeros(n_refine,1);



for i=1:size(Refinement,2)

   [err(i,:),n_dofs(i)]=IGA_3D_Poisson(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,Refinement(i),t, test_case);
 
 %  [err(i,:),n_dofs(i)]=IGA_3D_Poisson_Store_Basis(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,Refinement(i),t, test_case);

end



format short e
disp('The degree of NURBS basis is ')
disp(pu+t)

 disp('================================================================')
 disp(['The L2 norm error            ||       The H1 norm error'  ])
 
 disp('----------------------------------------------------------------')
 disp([err])
 disp('----------------------------------------------------------------')
disp('The #DOFs:')
disp(n_dofs)
disp('-----------------------------------------------------------------')
  disp('The convergence order in L2, H1 are:  ')
 
  disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2)])
 disp('================================================================')

end



