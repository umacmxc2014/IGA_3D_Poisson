function [err,u_h] = L2_project2dirichlet_bnd(ConPts,knotU,weights,p, Ubar,wbar, Ubreaks, uNoEs, u_d)

% 这个函数实现了把二阶PDE的非齐次狄利克雷边界函数 $u_d$ $L^2$ 投影到 边界上的 B-样条基函数张成的有限维空间.

wbar=wbar'; %注意到，这里的权系数需要是行向量.

addpath('./IGA_Grid_data/')
addpath('.//NURBS/')
addpath('./quadrature/')



ele_dof=p+1;

coordinate=zeros(uNoEs,2);
element=zeros(uNoEs,ele_dof);
knotSpanIndex=zeros(uNoEs,1);

for i=1:uNoEs
coordinate(i,:)=[Ubreaks(i),Ubreaks(i+1)];
span=findspan(Ubar,p,Ubreaks(i));
knotSpanIndex(i)=span;
element(i,:)=span-p:span;
end

dofs=length(wbar);
A=sparse(dofs,dofs);
rhs=zeros(dofs,1);


np=p+3 ;% 每一条边上的数值积分点的个数.


for i=1:uNoEs
	  row=element(i,:);
     
    ue=coordinate(i,:);
    span=knotSpanIndex(i);
    m_row = ele_dof;
    n_column = ele_dof + 1;
    x1={np,ue(1),ue(2),m_row, n_column  };
    s=Gauss_1d(@(u)quad_Ae_Fe(u,Ubar,wbar,ConPts,weights,knotU,p,u_d,span), x1{:});
     
    A(row,row)=A(row,row)+s(:,1:ele_dof);
    rhs(row)=rhs(row)+s(:,end);
end

u_h=A\rhs;

 err=0;

 for i=1:uNoEs
   row=element(i,:);
   ue=coordinate(i,:);
   span=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ue(1),ue(2),m_row,n_column};
   s=Gauss_1d(@(u)quad_err_L2(u,Ubar,wbar,ConPts,weights,knotU,p,u_d,u_h,span),x1{:});
   err=err+s;
   end
   err(1)=sqrt(err(1));

end


function s=quad_Ae_Fe(u,Ubar,wbar,ConPts,weights,U,p,u_d,span)

% 把非齐次狄利克雷边界上的边界函数 L^2 投影到边界上的NURBS空间。
%　这里计算L^2投影时，每一条边上所需要的质量矩阵和右端项。

[W,DW,C,DC]=NurbsCurve(ConPts,U,weights,p,u);
Uders=bspbasisDers(Ubar,p,u,0);
Nu = Uders(1,:);
Nu = Nu'; % 当前边上的 (p+1)个 非0的 B-样条基函数组成的列向量.

J=norm(DC,2);

Ae=Nu*Nu'*J; % Mass matrix.

Fe=u_d(C(1),C(2))*Nu*J;

s=[Ae,Fe];
end

function s=quad_err_L2(u,Ubar,wbar,ConPts,weights,U,p,u_exact,u_h,span)

[W,DW,C,DC]=NurbsCurve(ConPts,U,weights,p,u);
Uders=bspbasisDers(Ubar,p,u,0);
Nu=Uders(1,:);  % 当前边上的 (p+1)个 非0的 B-样条基函数组成的行向量.


J=norm(DC,2);

s=(u_exact(C(1),C(2)) - Nu*u_h(span-p:span)  ).^2*J;
% s_inner=(u_exact(C(1),C(2)) - Nu*u_h(span-p:span)  )*Nu*u_h(span-p:span)*J;
% s=[s,s_inner];

end
