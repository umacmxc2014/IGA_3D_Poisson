function [err,n_dofs]=IGA_3D_Poisson_Store_Basis(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,Refinement,t, test_case)

% 现在利用 NURBS 基函数 中的B-spline 基函数作为有限元空间。

%=====================
% Input:
% ConPts are the control points;
% weights are the weights;
% knotU contains the knot vector in the u direction;
% pu is the degree of NURBS surface in the u direction;
% knotV contains the knot vector in the v direction;
% pv is the degree of NURBS surface in the v direction;
% Refinement mens the times of the h-refinement of both u and v direction;
% t denotes the order to be elevated with analogous to p-refin, i.e., the
% ultimate degree of NURBS  basis functions is (pu+t) in u direction, and
% (pv+t) in the v direction;

tic

addpath('./IGA_Grid_data/')
addpath('./NURBS/')
addpath('./quadrature/')

 ConPts_0 =ConPts; 
 weights_0 = weights;
 knotU_0 = knotU;  
 knotV_0 = knotV;
 knotW_0 = knotW;
 pu_0 = pu;
 pv_0 = pv; 
 pw_0 = pw;



if(strcmp(test_case ,'rectangle'))

u_Exact=@(x,y,z)sin(pi*x)*sin(pi*y)*sin(pi*z);% Exact solution of Poisson equation;
f=@(x,y,z)  3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z); % The right hand side of Poisson equation;
u_d=@(x,y,z) 0;
u_Grad=@(x,y,z) pi*[cos(pi*x)*sin(pi*y)*sin(pi*z), sin(pi*x)*cos(pi*y)*sin(pi*z), sin(pi*x)*sin(pi*y)*cos(pi*z) ];

end

if(strcmp(test_case,'quarter'))
    
    
u_Exact=@(x,y)  x.^2.*y.^2*sin(pi*(x.^2+y.^2-2));% Exact solution of Poisson equation;
f=@(x,y) -(  20*pi*x.^2.*y.^2.*cos(pi*(x.^2+y.^2-2)) + ...
      2*(x^2+y^2-2*pi*pi*x^4*y^2-2*pi*pi*x^2*y^4)*sin(pi*(x.^2+y.^2-2))  );% The right hand side of Poisson equation;
u_d=@(x,y) x.^2.*y.^2*sin(pi*(x.^2+y.^2-2));

u_Grad=@(x,y) [ y^2*(2*x*sin(pi*(x.^2+y.^2-2)) +2*pi*x^3*cos(pi*(x.^2+y.^2-2))) , ...
               x^2*(2*y*sin(pi*(x.^2+y.^2-2)) +2*pi*y^3*cos(pi*(x.^2+y.^2-2))) ] ;
    

    
end


if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction, to (pv+t) in the v direction;
            % and to (pw+t) in the w-direction.
            
[Q,wbar,Ubar,Vbar,Wbar]=IGADegreeElevVolume(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,t);

 ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar; knotW=Wbar;
 pu=pu+t;  pv=pv+t; pw = pw+t; 
end

nurbsInfo=Iga_3d_grid(knotU,pu,knotV,pv,knotW,pw,weights,Refinement);








Element=nurbsInfo.Element;
Coordinate=nurbsInfo.Coordinate;

Ubar=nurbsInfo.Ubar;
Vbar=nurbsInfo.Vbar;
Wbar=nurbsInfo.Wbar;
N1=nurbsInfo.N1;
N2=nurbsInfo.N2;
N3=nurbsInfo.N3;

NoEs=nurbsInfo.NoEs;
n_dofs=nurbsInfo.n_dofs;



uNoEs=nurbsInfo.uNoEs;
vNoEs=nurbsInfo.vNoEs;
wNoEs=nurbsInfo.wNoEs;






[n_conpts_u,n_conpts_v,n_conpts_w,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.


% 我们现在只考虑 w=0 的边界上为非齐次狄利克雷边界，在其它边界上都是齐次狄利克雷边界条件.
UV_ConPts=zeros(n_conpts_u,n_conpts_v,DIM);

k=1;
for d=1:DIM
  for i=1:n_conpts_u
      for j=1:n_conpts_v
    UV_ConPts(i,j,d)=ConPts(i,j,k,d); % k=1对应于 w=0 的面.
end
  end
end


UV_weights = zeros(n_conpts_u,n_conpts_v);
  for i=1:n_conpts_u
      for j=1:n_conpts_v
          UV_weights(i,j) = weights(i,j,k);
      end
  end
  


 [err_bnd_L2_proj,u_d_h] = L2_3D_project2dirichlet_bnd(UV_ConPts,knotU,knotV,UV_weights,pu,pv,Ubar,Vbar,u_d);
 



% disp('The L2 projection error  is ')
% disp(u_d_h)




disp('The degree of the  NURBS basis is ')
disp(pu)
disp('The number of elements in the mesh and the number of DOFs are: ')
disp([NoEs,n_dofs])



Eledof=(pu+1)*(pv+1)*(pw+1); % The number of DOFs in an "Element";



 
nurbs_F = F_DF_Info(ConPts_0,weights_0,knotU_0,knotV_0,knotW_0,Ubar,pu_0,Vbar,pv_0,Wbar,pw_0,t);


nurbs_basis = Basis_Funcs_Info(Ubar,pu,Vbar,pv,Wbar,pw);

[A,rhs]=solve_laplace(nurbsInfo, nurbs_F, nurbs_basis,f);






% 我们现在考虑 w=0 (下底面对应的边界) 边界为 非齐次 Dirichlet 边界。

np = pu+1;

Element_w_0 = nurbsInfo.Element_w_0;

for e=1:(uNoEs*vNoEs) % Loop for the  elements lying on the boundary w=0;
    row=Element(e,:);
    ue=Coordinate(e,1:2);
    ve=Coordinate(e,3:4);
    we=Coordinate(e,5:6);
    uv_index=Element_w_0(e,:);
  
    m_row = Eledof;
    n_column = 1;

	x1={ConPts,weights,knotU,pu,knotV,pv,knotW,pw,Ubar,Vbar,Wbar,uv_index};
	x2={np,ue(1),ue(2),ve(1),ve(2),we(1),we(2), m_row, n_column };
    
    s=Gauss_3d(@(u,v,w)quad_Ae_Fe_uv(u,v,w,x1{:},u_d_h),x2{:});%----s=Gauss_2d(f,np,a,b,c,d)

    rhs(row)=rhs(row)-s;
end



[A,rhs]=Iga_3d_bc(A,rhs,N1,N2,N3);

% spy(A)

Uh=A\rhs;

% 这时候先把齐次解$u_0$算出来了。
% 现在把非齐次狄利克雷边界的解加上上面的非齐次解里。

Uh(1:N1*N2)=Uh(1:N1*N2)+u_d_h;


%%

% Compute the L2 error and the semi H1 error.
err = Err_L2(nurbsInfo, nurbs_F, nurbs_basis,Uh,u_Exact, u_Grad);

%%
disp('The total simulation time is')
toc

disp('===========================================')
    
end



function [A,rhs]=solve_laplace(nurbsInfo, nurbs_F, nurbs_basis,f)


n_dofs = nurbsInfo.n_dofs;
% A = sparse(n_dofs,n_dofs);
rhs =zeros(n_dofs,1);
NoEs = nurbsInfo.NoEs;
Element=nurbsInfo.Element;
F   = nurbs_F.ele_F;
DF = nurbs_F.ele_DF;  
Jacobian = nurbs_F.Jacobian;

basis_funcs = nurbs_basis.ele_basis_funcs;          
basis_grad  = nurbs_basis.ele_basis_grad;

n_gps = nurbs_basis.n_gps;
n_ele_dofs =  nurbs_basis.n_ele_dofs;
Ae = zeros(n_ele_dofs,n_ele_dofs);
Fe = zeros(n_ele_dofs,1);

row_index = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
column_index = row_index;
value = row_index;
global_index = 1;

for e = 1:NoEs
    row = Element(e,:);
    Ae = 0*Ae; Fe = 0*Fe;
    ele_F = F{e};
    ele_basis_funcs = basis_funcs{e};
    for i=1:n_gps
        i_ele_basis_funcs = ele_basis_funcs(:,i);
        i_F=ele_F(:,i);
        ele_DF =DF{e,i};
        ele_basis_grad = basis_grad{e,i}/ele_DF;
        ele_Jacobi = Jacobian(e,i);
        Ae = Ae + ele_Jacobi*ele_basis_grad*(ele_basis_grad'); 
        Fe = Fe + f(i_F(1),i_F(2),i_F(3))*i_ele_basis_funcs *ele_Jacobi;
    end
       % A(row,row) = A(row,row) + Ae;
       rhs(row) = rhs(row) + Fe;
       for i1=1:n_ele_dofs
           for j1=1:n_ele_dofs
               row_index(global_index) = row(i1);
               column_index(global_index) = row(j1);
               value(global_index) =  value(global_index) + Ae(i1,j1);
               global_index = global_index + 1;
           end
       end
end

A =sparse(row_index,column_index,value,n_dofs,n_dofs); 

end




function s=quad_Ae_Fe_uv(u,v,w,ConPts,weights,knotU,pu,knotV,pv,knotW,pw,Ubar,Vbar,Wbar,uv_index,u_d_h )
% 这个函数是为了计算 边界函数 $u_d$ 的 $L^2$ 投影的梯度与 $v_h$的梯度的积分。
% ----- Eledof=(pu+1)*(pv+1)*(pw+1);
[F,DF]=NurbsVolume(ConPts,weights,knotU,pu,u,knotV,pv,v,knotW,pw,w);

Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pv,v,1);
Wders=bspbasisDers(Wbar,pw,w,1);
Nu=Uders(1,:)';    DNu=Uders(2,:)';
Nv=Vders(1,:);     DNv=Vders(2,:);
Nw=Wders(1,:);   DNw=Wders(2,:);
J=abs(det(DF));


n_ele_dofs_uv = (pu+1)*(pv+1);

B_uv = Nu*Nv;
B_uv = B_uv(:);

DBu_v = DNu*Nv;
DBu_v = DBu_v(:);
DBu = DBu_v*Nw;
DBu = DBu(:);
% This array stores the first component of the gradient of all basis functions B_{ijk}(u,v,w),  i.e., $\partial B_{ijk}(u,v,w)/\partial u$,
% where $B_{ijk}=N_{i,pu}(u)M_{j,pv}(v)L_{k,pw}(w)$.


u_DBv = Nu*DNv;
u_DBv = u_DBv(:);
DBv = u_DBv*Nw;
DBv = DBv(:);
% This array stores the second component of the gradient of all basis functions B_{ijk}(u,v,w),  i.e., $\partial B_{ijk}(u,v,w)/\partial v$,
% where $B_{ijk}=N_{i,pu}(u)M_{j,pv}(v)L_{k,pw}(w)$.



DBw = B_uv*DNw;
DBw = DBw(:);
% This array stores the third component of the gradient of all basis functions B_{ijk}(u,v,w),  i.e., $\partial B_{ijk}(u,v,w)/\partial w$,
% where $B_{ijk}=N_{i,pu}(u)M_{j,pv}(v)L_{k,pw}(w)$.


DB = [DBu,DBv,DBw]/DF; % 当前非 0 的 B-样条基函数的梯度组成的矩阵.

DB_w_0 = DB(1:n_ele_dofs_uv,:);  % $N_{i,pu}(u)M_{j,pv}(v)L_{1,pw}(w)$.

D_u_d_h_w0= u_d_h(uv_index)'*DB_w_0*J; D_u_d_h_w0 = D_u_d_h_w0'; % 这时候是列向量了.
s = DB*D_u_d_h_w0;




end


function s=Err_L2(nurbsInfo, nurbs_F, nurbs_basis,Uh,u_Exact, u_Grad)

 
NoEs = nurbsInfo.NoEs;
Element=nurbsInfo.Element;
F   = nurbs_F.ele_F;
DF = nurbs_F.ele_DF;  
Jacobian = nurbs_F.Jacobian;

basis_funcs = nurbs_basis.ele_basis_funcs;          
basis_grad  = nurbs_basis.ele_basis_grad;

n_gps = nurbs_basis.n_gps;


err_L2 = 0;
err_semi_H1 = 0;

for e = 1:NoEs
    row = Element(e,:);
    ele_F = F{e};
    ele_basis_funcs = basis_funcs{e};
    for i=1:n_gps
        i_ele_basis_funcs = ele_basis_funcs(:,i);
        i_F=ele_F(:,i);
        ele_DF =DF{e,i};
        ele_basis_grad = basis_grad{e,i}/ele_DF;
        ele_Jacobi = Jacobian(e,i);
        err_L2 = err_L2 + (u_Exact(i_F(1),i_F(2),i_F(3)) -  Uh(row)'*i_ele_basis_funcs)^2*ele_Jacobi;
        err_semi_H1 = err_semi_H1 + sum((u_Grad(i_F(1),i_F(2),i_F(3)) - Uh(row)'*ele_basis_grad ).^2)*ele_Jacobi;
    end

end


s=sqrt([err_L2,err_semi_H1]);


end


