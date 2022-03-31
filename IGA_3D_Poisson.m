function [err,n_dofs]=IGA_3D_Poisson(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,Refinement,t, test_case)

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

nurbs_original.ConPts = ConPts;
nurbs_original.weights = weights;
nurbs_original.knotU = knotU; 
nurbs_original.knotV  = knotV; 
nurbs_original.knotW = knotW; 
nurbs_original.pu = pu;
nurbs_original.pv = pv;
nurbs_original.pw = pw;


addpath('./IGA_Grid_data/')
addpath('./NURBS/')
addpath('./quadrature/')



if(strcmp(test_case ,'rectangle'))

u_Exact=@(x,y,z)sin(pi*x)*sin(pi*y)*sin(pi*z);% Exact solution of Poisson equation;
f=@(x,y,z)  3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z); % The right hand side of Poisson equation;
u_d=@(x,y,z) 0;
u_Grad=@(x,y,z) pi*[cos(pi*x)*sin(pi*y)*sin(pi*z), sin(pi*x)*cos(pi*y)*sin(pi*z), sin(pi*x)*sin(pi*y)*cos(pi*z) ];

end

if(strcmp(test_case,'cylinder'))
    
u_Exact=@(x,y,z)sin(pi*(x^2+y^2-1))*sin(z-1);% Exact solution of Poisson equation;
f=@(x,y,z) (4*pi*pi*x^2+4*pi*pi*y^2+1)*sin(pi*(x^2+y^2-1))*sin(z-1) -4*pi*cos(pi*(x^2+y^2-1))*sin(z-1) ; % The right hand side of Poisson equation;
u_d=@(x,y,z) sin(pi*(x^2+y^2-1))*sin(z-1);

u_Grad=@(x,y,z) [2*pi*x*cos(pi*(x^2+y^2-1))*sin(z-1),2*pi*y*cos(pi*(x^2+y^2-1))*sin(z-1),...
                             sin(pi*(x^2+y^2-1))*cos(z-1) ] ;
    

    
end


if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction, to (pv+t) in the v direction;
            % and to (pw+t) in the w-direction.
            
[Q,wbar,Ubar,Vbar,Wbar]=IGADegreeElevVolume(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,t);

 ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar; knotW=Wbar;
 pu=pu+t;  pv=pv+t; pw = pw+t; 
end

nurbs_refine=Iga_3d_grid(knotU,pu,knotV,pv,knotW,pw,weights,Refinement);


Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
Wbar=nurbs_refine.Wbar;
N1=nurbs_refine.N1;
N2=nurbs_refine.N2;
N3=nurbs_refine.N3;

NoEs=nurbs_refine.NoEs;
n_dofs=nurbs_refine.n_dofs;




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
  


 [err,u_d_h] = L2_3D_project2dirichlet_bnd(UV_ConPts,knotU,knotV,UV_weights,pu,pv,Ubar,Vbar,u_d);
 

 disp('L2 projection error is ')
 disp(err)


% disp('The L2 projection error  is ')
% disp(u_d_h)



% if pu==1
% np=pu ;% The number of Gauss quadrature points in  element;
% else
    np = pu + 1 ;
% end

if(np>=9)
    np=9;
end


disp('The degree of the  NURBS basis is ')
disp(pu)
disp('The number of elements in the mesh and the number of DOFs are: ')
disp([NoEs,n_dofs])



% disp('*************************************')
% disp('The time for assembing the stiffness matrix and right hand side:')
% tic
[A,rhs] = solve_laplace_A_F(nurbs_original,nurbs_refine,f);
% toc
% disp('*************************************')


% 我们现在考虑 w=0 (下底面对应的边界) 边界为 非齐次 Dirichlet 边界。


 rhs = modify_rhs_w_0(nurbs_original,nurbs_refine,u_d_h,rhs);





[A,rhs]=Iga_3d_bc(A,rhs,N1,N2,N3);

% spy(A)

Uh=A\rhs;

% 这时候已经把齐次解$u_0$算出来了。
% 下面把非齐次狄利克雷边界的解加上上面的非齐次解里。

Uh(1:N1*N2)=Uh(1:N1*N2)+u_d_h;

% disp('---------------------------------------------------------------------------')
% disp('The time for computing the L2 and semi-H1 norm error is')
% tic
err=Compute_Error(nurbs_original,nurbs_refine,Uh,u_Exact, u_Grad);
% toc
% disp('---------------------------------------------------------------------------')
disp('The total simulation time is')
toc

disp('===========================================')

end




