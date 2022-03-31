function rhs = modify_rhs_w_0(nurbs_original,nurbs_refine,u_d_h,rhs)


u_d_h = u_d_h';

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
knotW_o   = nurbs_original.knotW; 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;
pw_o        =  nurbs_original.pw;




Element=nurbs_refine.Element;
knotU=nurbs_refine.Ubar;
knotV=nurbs_refine.Vbar;
knotW=nurbs_refine.Wbar;
UBreaks=nurbs_refine.UBreaks;   % u 方向上节点向量中的断点.
VBreaks=nurbs_refine.VBreaks;    % v 方向上节点向量中的断点.
WBreaks=nurbs_refine.WBreaks;    % w 方向上节点向量中的断点.



uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;
% wNoEs=nurbs_refine.wNoEs;

pu            =  nurbs_refine.pu;
pv            =  nurbs_refine.pv;
pw            =  nurbs_refine.pw;

Element_w_0 = nurbs_refine.Element_w_0;


u_np = pu + 1;   % It seems that if pu >=3, u_np = pu is enough!!!
v_np = pv + 1;   
w_np = pw+1; 

[gp_u,gw_u]=quad_info_1d(u_np); % The quadrature points and quadrature weights in the [-1,1].

u_ele_basis_funcs_o =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);
u_ele_basis_grad_o  =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);

u_knot_span_o = zeros(uNoEs,1);
u_jacobian    = zeros(uNoEs, u_np);



basis = zeros(pu_o+1,u_np);
basis_grad = basis;
for i=1:uNoEs
	ue = UBreaks(i:i+1);
    uJ=(ue(2)-ue(1))/2;
    u_knot_span_o(i) = findspan(knotU_o,pu_o,ue(1));
    for j=1:u_np
     u = uJ*gp_u(j)+(ue(1)+ue(2))/2;
     Uders=bspbasisDers(knotU_o,pu_o,u,1); 
     basis(:,j) = Uders(1,:);
     basis_grad(:,j) = Uders(2,:);
     u_jacobian(i,j) = uJ*gw_u(j);
    end
u_ele_basis_funcs_o{i}=basis;
u_ele_basis_grad_o{i}=basis_grad;
end

u_ele_basis_funcs =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);
u_ele_basis_grad  =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);
u_knot_span = zeros(uNoEs,1);
basis = zeros(pu+1,u_np);
basis_grad = basis;
for i=1:uNoEs
	ue = UBreaks(i:i+1);
    uJ=(ue(2)-ue(1))/2;
    u_knot_span(i) = findspan(knotU,pu,ue(1));
    for j=1:u_np
     u = uJ*gp_u(j)+(ue(1)+ue(2))/2;
     Uders=bspbasisDers(knotU,pu,u,1); 
     basis(:,j) = Uders(1,:);
     basis_grad(:,j) = Uders(2,:);
    end
u_ele_basis_funcs{i}=basis;
u_ele_basis_grad{i}=basis_grad;
end



[gp_v,gw_v]=quad_info_1d(v_np); % The quadrature points and quadrature weights in the [-1,1].

v_ele_basis_funcs_o =  cell(vNoEs,1); % zeros(vNoEs,pv+1,v_np);
v_ele_basis_grad_o  =  cell(vNoEs,1); % zeros(vNoEs,pv+1,v_np);
v_knot_span_o = zeros(vNoEs,1);
v_jacobian    = zeros(vNoEs, v_np);

basis = zeros(pv_o+1,v_np);
basis_grad = basis;
for i=1:vNoEs
	ve = VBreaks(i:i+1);
    vJ=(ve(2)-ve(1))/2;
    v_knot_span_o(i) = findspan(knotV_o,pv_o,ve(1));
    for j=1:v_np
     v = vJ*gp_v(j)+(ve(1)+ve(2))/2;
     Vders=bspbasisDers(knotV_o,pv_o,v,1); 
     basis(:,j) = Vders(1,:);
     basis_grad(:,j) = Vders(2,:);
     v_jacobian(i,j) = vJ*gw_v(j);
    end
v_ele_basis_funcs_o{i}=basis;
v_ele_basis_grad_o{i}=basis_grad;
end

basis = zeros(pv+1,v_np);
basis_grad = basis;
v_ele_basis_funcs =  cell(vNoEs,1); % zeros(uNoEs,pu+1,u_np);
v_ele_basis_grad  =  cell(vNoEs,1); % zeros(uNoEs,pu+1,u_np);
v_knot_span = zeros(vNoEs,1);
for i=1:vNoEs
	ve = VBreaks(i:i+1);
    vJ=(ve(2)-ve(1))/2;
    v_knot_span(i) = findspan(knotV,pv,ve(1));
    for j=1:v_np
     v = vJ*gp_v(j)+(ve(1)+ve(2))/2;
     Vders=bspbasisDers(knotV,pv,v,1); 
     basis(:,j) = Vders(1,:);
     basis_grad(:,j) = Vders(2,:);
    end
v_ele_basis_funcs{i}=basis;
v_ele_basis_grad{i}=basis_grad;
end



% 由于我们现在只考虑在w=0的边界上为非齐次狄利边界条件，所以我们只需要考虑[w_1,w_{1+p+1}]上的非０基函数.

[gp_w,gw_w]=quad_info_1d(w_np); % The quadrature points and quadrature weights in the [-1,1].

w_ele_basis_funcs_o =  zeros(pw_o+1,w_np);  
w_ele_basis_grad_o  =  zeros(pw_o+1,w_np);  

w_jacobian    = zeros(1, w_np);


    i=1;
	we = WBreaks(i:i+1);
    wJ=(we(2)-we(1))/2;
    w_knot_span_o = findspan(knotW_o,pw_o,we(1));
    for j=1:w_np
     w = wJ*gp_w(j)+(we(1)+we(2))/2;
     Wders=bspbasisDers(knotW_o,pw_o,w,1); 
     w_ele_basis_funcs_o(:,j) = Wders(1,:);
     w_ele_basis_grad_o(:,j) = Wders(2,:);
     w_jacobian(j) = wJ*gw_w(j);
    end




w_ele_basis_funcs =  zeros(pw+1,w_np);
w_ele_basis_grad  =  zeros(pw+1,w_np);

     i=1;
	we = WBreaks(i:i+1);
    wJ=(we(2)-we(1))/2;
    for j=1:w_np
     w = wJ*gp_w(j)+(we(1)+we(2))/2;
     Wders=bspbasisDers(knotW,pw,w,1); 
     w_ele_basis_funcs(:,j) = Wders(1,:);
     w_ele_basis_grad(:,j) = Wders(2,:);
    end





n_gps = u_np*v_np *w_np;


n_ele_dofs = (pu+1)*(pv+1)*(pw+1);

n_ele_dofs_w_0 = (pu+1)*(pv+1); % The number of DOFs in each element of the face w=0

Fe = zeros(n_ele_dofs,1);


DF = cell(n_gps,1);
Jacobian = zeros(n_gps,1);





    for j=1:vNoEs
        for i=1:uNoEs
   
            e = i + (j-1)*uNoEs;
            row = Element(e,:);
            
            Fe = 0*Fe;
            
            uv_index = Element_w_0(e,:);

            u_basis_i_o   =  u_ele_basis_funcs_o{i};
            v_basis_j_o   =  v_ele_basis_funcs_o{j};
            w_basis_k_o =  w_ele_basis_funcs_o;
            uv_basis_ij_o = kron(v_basis_j_o,u_basis_i_o);
            uvw_basis_ijk_o = kron(w_basis_k_o,uv_basis_ij_o);
            
            Du_basis_i_o =  u_ele_basis_grad_o{i};
            Du_v_basis_ij_o = kron(v_basis_j_o,Du_basis_i_o);
            Du_vw_basis_ijk_o = kron(w_basis_k_o,Du_v_basis_ij_o);
            
            Dv_basis_j_o =  v_ele_basis_grad_o{j};
            Dv_u_basis_ij_o = kron(Dv_basis_j_o,u_basis_i_o);
            Dv_uw_basis_ijk_o = kron(w_basis_k_o,Dv_u_basis_ij_o);  
            
            Dw_basis_k_o =  w_ele_basis_grad_o;
            Dw_uv_basis_ijk_o = kron(Dw_basis_k_o,uv_basis_ij_o);
            
            
            uspan_o  =  u_knot_span_o(i);
            vspan_o  =  v_knot_span_o(j);
            wspan_o =  w_knot_span_o;
            
            w = weights_o(uspan_o-pu_o:uspan_o,vspan_o-pv_o:vspan_o,wspan_o-pw_o:wspan_o);
            w = w(:);

            w_funcs_o = w'*uvw_basis_ijk_o; % 现在算出了权重函数在 u_np*v_np*w_np 个积分点处的值
            
            
            P_x = ConPts_o(uspan_o-pu_o:uspan_o,vspan_o-pv_o:vspan_o,wspan_o-pw_o:wspan_o,1);
            P_x = P_x(:);
            w_P_x = (w.*P_x)';
            F_x = (w_P_x* uvw_basis_ijk_o)./w_funcs_o;
            
            
            P_y = ConPts_o(uspan_o-pu_o:uspan_o,vspan_o-pv_o:vspan_o,wspan_o-pw_o:wspan_o,2);
            P_y = P_y(:);
            w_P_y = (w.*P_y)';
            F_y = (w_P_y* uvw_basis_ijk_o)./w_funcs_o;
            
            P_z = ConPts_o(uspan_o-pu_o:uspan_o,vspan_o-pv_o:vspan_o,wspan_o-pw_o:wspan_o,3);
            P_z = P_z(:);
            w_P_z = (w.*P_z)';
            F_z = (w_P_z* uvw_basis_ijk_o)./w_funcs_o;
            
            
            
            
            Dw_u_o   = w'*Du_vw_basis_ijk_o; 
            Dw_v_o   = w'*Dv_uw_basis_ijk_o; 
            Dw_w_o  = w'*Dw_uv_basis_ijk_o; 
            
        
            DF_x_u  = (w_P_x*Du_vw_basis_ijk_o - F_x.*Dw_u_o )./w_funcs_o;
            DF_x_v  = (w_P_x*Dv_uw_basis_ijk_o - F_x.*Dw_v_o )./w_funcs_o;
            DF_x_w = (w_P_x*Dw_uv_basis_ijk_o - F_x.*Dw_w_o )./w_funcs_o;

            DF_y_u  = (w_P_y*Du_vw_basis_ijk_o - F_y.*Dw_u_o )./w_funcs_o;
            DF_y_v  = (w_P_y*Dv_uw_basis_ijk_o - F_y.*Dw_v_o )./w_funcs_o;
            DF_y_w = (w_P_y*Dw_uv_basis_ijk_o - F_y.*Dw_w_o )./w_funcs_o;
            
            DF_z_u  = (w_P_z*Du_vw_basis_ijk_o - F_z.*Dw_u_o )./w_funcs_o;
            DF_z_v  = (w_P_z*Dv_uw_basis_ijk_o - F_z.*Dw_v_o )./w_funcs_o;
            DF_z_w = (w_P_z*Dw_uv_basis_ijk_o - F_z.*Dw_w_o )./w_funcs_o;
            
            
            u_gw = u_jacobian(i,:);         u_gw = u_gw(:);
            v_gw = v_jacobian(j,:);         v_gw  = v_gw(:);
            w_gw = w_jacobian;            w_gw = w_gw(:);
            uv_gw = kron(v_gw,u_gw); 
            uvw_gw = kron(w_gw,uv_gw);
            
            for k1=1:n_gps
            DF{k1}=[DF_x_u(k1)  DF_x_v(k1)  DF_x_w(k1);... 
                            DF_y_u(k1)  DF_y_v(k1)  DF_y_w(k1);...
                            DF_z_u(k1)  DF_z_v(k1)  DF_z_w(k1)];
            Jacobian(k1) = abs(det(DF{k1}))*uvw_gw(k1);
            end
            

            
            
            
            
            u_basis_i   =  u_ele_basis_funcs{i};
            v_basis_j   =  v_ele_basis_funcs{j};
            w_basis_k =  w_ele_basis_funcs;
            uv_basis_ij = kron(v_basis_j,u_basis_i);
           
            
            Du_basis_i  = u_ele_basis_grad{i};
            Du_v_basis_ij = kron(v_basis_j,Du_basis_i);
            Du_vw_basis_ijk = kron(w_basis_k,Du_v_basis_ij);
            
            Dv_basis_j  = v_ele_basis_grad{j};
            Dv_u_basis_ij = kron(Dv_basis_j,u_basis_i);
            Dv_uw_basis_ijk = kron(w_basis_k,Dv_u_basis_ij);
            
            Dw_basis_k  = w_ele_basis_grad;
            Dw_uv_basis_ijk = kron(Dw_basis_k,uv_basis_ij);
            

            
            for k1=1:n_gps
                basis_grad = [Du_vw_basis_ijk(:,k1),Dv_uw_basis_ijk(:,k1),Dw_uv_basis_ijk(:,k1)]/DF{k1};
                basis_grad_w_0 = basis_grad(1:n_ele_dofs_w_0,:); % The grad of all basis functions on the face w=0;
                grad_u_d_h_w_0 = u_d_h(uv_index)*basis_grad_w_0;
                Fe = Fe + basis_grad*(grad_u_d_h_w_0')*Jacobian(k1);
            end
            

            rhs(row) = rhs(row) - Fe;
                                           

        end
    end
    

end