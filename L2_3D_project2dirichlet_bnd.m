function [err, u_h] = L2_3D_project2dirichlet_bnd(ConPts,knotU,knotV,weights,pu,pv,Ubar,Vbar, u_d)

addpath('./IGA_Grid_data/')
addpath('.//NURBS/')
addpath('./quadrature/')


DIM=2;

% 这个函数实现了把二阶PDE的非齐次狄利克雷边界函数 $u_d$ $L^2$ 投影到 边界上的 B-样条基函数张成的有限维空间.

N1 = length(Ubar) - pu -1;   % The number of basis functions in the u-direction.
N2 = length(Vbar) - pv -1;   %  The number of basis functions in the v-direction.

uBreaks = unique(Ubar);  uNoEs = length(uBreaks) - 1;
vBreaks = unique(Vbar);  vNoEs =  length(vBreaks) - 1;

n_ele_dof=(pu+1)*(pv+1);

NoEs = uNoEs*vNoEs;

coordinate=zeros(NoEs,2*DIM);
ele_dof=zeros(NoEs,n_ele_dof);
knotSpanIndex=zeros(NoEs,DIM);

ind=1;
for j=1:vNoEs 
    for i=1:uNoEs
coordinate(ind,:)=[uBreaks(i),uBreaks(i+1), vBreaks(j),vBreaks(j+1) ];
u_span=findspan(Ubar,pu,uBreaks(i));  v_span=findspan(Vbar,pv,vBreaks(j));
knotSpanIndex(ind,:)=[u_span,v_span];

           local_index = 1;
   for j1=(v_span - pv):v_span
       for i1=(u_span - pu):u_span
           global_index = i1 + (j1-1)*N1;
           ele_dof(ind,local_index)=global_index;
           local_index = local_index + 1;
       end
   end
   ind = ind +1;
    end
end

n_dofs=N1*N2;
% A=sparse(n_dofs,n_dofs);
rhs=zeros(n_dofs,1);


np=pu+1;  % 每一条边上的数值积分点的个数.

row_index = zeros(NoEs*n_ele_dof*n_ele_dof,1);
column_index =row_index;
value = row_index;
global_index = 1;

for e = 1:NoEs
	row=ele_dof(e,:);
    ue=coordinate(e,:);
    m_row = n_ele_dof;
    n_column = n_ele_dof + 1;
    x1={np,ue(1),ue(2),ue(3),ue(4), m_row, n_column};
    s=Gauss_2d(@(u,v)quad_Ae_Fe(u,v,ConPts,weights,knotU,knotV,Ubar,Vbar,pu,pv,u_d), x1{:});
    Ae = s(:,1:n_ele_dof);
   % A(row,row)=A(row,row)+s(:,1:n_ele_dof);
   for i1=1:n_ele_dof
       for j1=1:n_ele_dof
           value(global_index) = value(global_index) + Ae(i1,j1);
           row_index(global_index) = row(i1);
           column_index(global_index) = row(j1);
           global_index = global_index + 1;
       end
   end
   
   
   
    rhs(row)=rhs(row)+s(:,end);
end

A =sparse(row_index,column_index,value,n_dofs,n_dofs);  

u_h=A\rhs;

  err=0;

  for e=1:NoEs
   row=ele_dof(e,:);
   ue=coordinate(e,:);
   m_row = 1;
   n_column = 1;
   x1={np,ue(1),ue(2),ue(3),ue(4), m_row,n_column};
   s=Gauss_2d(@(u,v)quad_err_L2(u,v,ConPts,weights,knotU,knotV, pu,pv, Ubar,Vbar,row,u_d,u_h),x1{:});
  err=err+s;
  end
   err(1)=sqrt(err(1));

end


function s=quad_Ae_Fe(u,v,ConPts,weights,knotU,knotV,Ubar,Vbar,pu,pv,u_d)

% 把非齐次狄利克雷边界上的边界函数 L^2 投影到边界上的NURBS空间。
%　这里计算L^2投影时，每一条边上所需要的质量矩阵和右端项。

[S,DF]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);

Uders=bspbasisDers(Ubar,pu,u,0);
Nu = Uders(1,:);
Nu = Nu'; %在点u处的 (pu+1)个 非0的 B-样条基函数组成的列向量.

Vders=bspbasisDers(Vbar,pv,v,0);
Nv = Vders(1,:); %在点v处的 (pv+1)个 非0的 B-样条基函数组成的列向量.

basis_funcs = Nu*Nv;
basis_funcs = basis_funcs(:);

A = DF(2,1)*DF(3,2) - DF(2,2)*DF(3,1);
B = DF(2,1)*DF(1,2) - DF(1,1)*DF(3,2);
C = DF(1,1)*DF(2,2) - DF(1,2)*DF(2,1);

J=sqrt(A*A+B*B+C*C);

Ae=basis_funcs*basis_funcs'*J; % Mass matrix.

Fe=u_d(S(1),S(2),S(3))*basis_funcs*J;

s=[Ae,Fe];
end

function s=quad_err_L2(u,v,ConPts,weights,knotU,knotV, pu,pv, Ubar,Vbar,row,u_d,u_h)

[S,DF]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);

Uders=bspbasisDers(Ubar,pu,u,0);
Nu=Uders(1,:);  % 当前边上的 (p+1)个 非0的 B-样条基函数组成的行向量.
Nu=Nu';

Vders=bspbasisDers(Vbar,pv,v,0);
Nv = Vders(1,:); %在点v处的 (pv+1)个 非0的 B-样条基函数组成的列向量.

basis_funcs = Nu*Nv;
basis_funcs = basis_funcs(:);
basis_funcs = basis_funcs'; % 化为行向量

A = DF(2,1)*DF(3,2) - DF(2,2)*DF(3,1);
B = DF(2,1)*DF(1,2) - DF(1,1)*DF(3,2);
C = DF(1,1)*DF(2,2) - DF(1,2)*DF(2,1);

J=sqrt(A*A+B*B+C*C);

s=(u_d(S(1),S(2),S(3)) - basis_funcs *u_h(row)  ).^2*J;


end
