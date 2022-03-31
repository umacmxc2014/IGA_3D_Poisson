function nurbsInfo = Basis_Funcs_Info(Ubar,pu,Vbar,pv,Wbar,pw)

UBreaks=unique(Ubar);   % u 方向上节点向量中的断点.
VBreaks=unique(Vbar);   % v 方向上节点向量中的断点.
WBreaks=unique(Wbar);   % w 方向上节点向量中的断点.

uNoEs=length(UBreaks)-1;      %  u 方向上的区间数。
vNoEs=length(VBreaks)-1;      %  v 方向上的区间数。
wNoEs=length(WBreaks)-1;      %  w 方向上的区间数。
NoEs=uNoEs*vNoEs*wNoEs;       %  计算区域上的区间总数。

 

DIM=3;
n_ele_dofs = (pu+1)*(pv+1)*(pw+1);


% tic

u_np = pu + 1;              % The number of Gauss quadrature points in the u-direction
[gp,gw]=quad_info_1d(u_np); % The quadrature points and quadrature weights in the [-1,1].


u_ele_basis_funcs = cell(uNoEs,1);    % zeros(uNoEs,pu+1,u_np);
u_ele_basis_grad  = cell(uNoEs,1);    % zeros(uNoEs,pu+1,u_np);
basis = zeros(pu+1,u_np);
basis_grad = basis;


for i=1:uNoEs
    ue = UBreaks(i:i+1);
    uJ=(ue(2)-ue(1))/2;
    for j=1:u_np
     u = uJ*gp(j)+(ue(1)+ue(2))/2;
     Uders=bspbasisDers(Ubar,pu,u,1); 
     basis(:,j) = Uders(1,:);
     basis_grad(:,j) = Uders(2,:);
    end
u_ele_basis_funcs{i} = basis;
u_ele_basis_grad{i} = basis_grad;
end

v_np = pv + 1;
[gp,gw]=quad_info_1d(v_np); % The quadrature points and quadrature weights in the [-1,1].
v_ele_basis_funcs =   cell(vNoEs,1); %  zeros(vNoEs,pv+1,v_np);
v_ele_basis_grad  =   cell(vNoEs,1); %  zeros(vNoEs,pv+1,v_np);
basis = zeros(pv+1,v_np);
basis_grad = basis;


for i=1:vNoEs
	ve = VBreaks(i:i+1);
    vJ=(ve(2)-ve(1))/2;
    for j=1:v_np
     v = vJ*gp(j)+(ve(1)+ve(2))/2;
     Vders=bspbasisDers(Vbar,pv,v,1); 
     basis(:,j) = Vders(1,:);
     basis_grad(:,j) = Vders(2,:);
    end
    v_ele_basis_funcs{i} = basis;
    v_ele_basis_grad{i} = basis_grad;
end

w_np = pw + 1;
[gp,gw]=quad_info_1d(w_np); % The quadrature points and quadrature weights in the [-1,1].
w_ele_basis_funcs = cell(wNoEs,1);  % zeros(wNoEs,pw+1,w_np);
w_ele_basis_grad  = cell(wNoEs,1);  % zeros(wNoEs,pw+1,w_np);
basis = zeros(pw+1,w_np);
basis_grad = basis;

for i=1:wNoEs
	we = WBreaks(i:i+1);
    wJ=(we(2)-we(1))/2;
    for j=1:w_np
     w = wJ*gp(j)+(we(1)+we(2))/2;
     Wders=bspbasisDers(Wbar,pw,w,1); 
     basis(:,j) = Wders(1,:);
     basis_grad(:,j) = Wders(2,:);
    end
    w_ele_basis_funcs{i} = basis;
    w_ele_basis_grad{i} = basis_grad;
end

n_gps = u_np*v_np*w_np;


ele_basis_funcs =  cell(NoEs,1);        %  zeros(NoEs,n_ele_dofs,n_gps);
ele_basis_grad  =  cell(NoEs,n_gps);  %  zeros(NoEs,n_gps,n_ele_dofs,DIM);
 

for k=1:wNoEs
	for j=1:vNoEs
		for i=1:uNoEs
			e = i + (j-1)*uNoEs + (k-1)*uNoEs*vNoEs;
            u_basis_i =  u_ele_basis_funcs{i};
            v_basis_j =  v_ele_basis_funcs{j};
            w_basis_k = w_ele_basis_funcs{k};
            uv_basis_ij = kron(v_basis_j,u_basis_i);
            uvw_basis_ijk = kron(w_basis_k,uv_basis_ij);
            ele_basis_funcs{e} = uvw_basis_ijk;
            
            Du_basis_i =  u_ele_basis_grad{i};
            Du_v_basis_ij = kron(v_basis_j,Du_basis_i);
            Du_vw_basis_ijk = kron(w_basis_k,Du_v_basis_ij);
             
            
            Dv_basis_j =  v_ele_basis_grad{j};
            Dv_u_basis_ij = kron(Dv_basis_j,u_basis_i);
            Dv_uw_basis_ijk = kron(w_basis_k,Dv_u_basis_ij);
             
            
            
            Dw_basis_k =  w_ele_basis_grad{k};
            Dw_uv_basis_ijk = kron(Dw_basis_k,uv_basis_ij);
           
            
            for k1=1:n_gps
                ele_basis_grad{e,k1} = [Du_vw_basis_ijk(:,k1), Dv_uw_basis_ijk(:,k1), Dw_uv_basis_ijk(:,k1)];
            end
            
        end
    end
end








 nurbsInfo.ele_basis_funcs = ele_basis_funcs;          
 nurbsInfo.ele_basis_grad = ele_basis_grad;
 nurbsInfo.n_gps = n_gps;
 nurbsInfo.n_ele_dofs = n_ele_dofs;


end

