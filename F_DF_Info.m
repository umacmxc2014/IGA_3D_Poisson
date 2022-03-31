function nurbsInfo = F_DF_Info(ConPts,weights,knotU,knotV,knotW,Ubar,pu,Vbar,pv,Wbar,pw,t)


UBreaks=unique(Ubar);       % 当前计算网格下u 方向上节点向量中的断点
VBreaks=unique(Vbar);       % 当前计算网格下v 方向上节点向量中的断点
WBreaks=unique(Wbar);     % 当前计算网格下w 方向上节点向量中的断点

uNoEs=length(UBreaks)-1;          %   当前计算网格下u 方向上的区间数
vNoEs=length(VBreaks)-1;          %   当前计算网格下v 方向上的区间数
wNoEs=length(WBreaks)-1;        %   当前计算网格下w 方向上的区间数
NoEs=uNoEs*vNoEs*wNoEs;       %   当前计算网格下计算区域上的区间总数



DIM=3;





u_np = pu + t + 1; % The number of of Gauss quadrature points in the u-direction
[gp,gw]=quad_info_1d(u_np); % The quadrature points and quadrature weights in the [-1,1].


u_ele_basis_funcs = cell(uNoEs,1);  %  zeros(uNoEs,pu+1,u_np);
u_ele_basis_grad  =  cell(uNoEs,1); %  zeros(uNoEs,pu+1,u_np);
u_knot_span = zeros(uNoEs,1);
u_jacobian = zeros(uNoEs, u_np);
basis = zeros(pu+1,u_np);
basis_grad = basis;

for i=1:uNoEs
	ue = UBreaks(i:i+1);
    uJ=(ue(2)-ue(1))/2;
    u_knot_span(i) = findspan(knotU,pu,UBreaks(1));
    for j=1:u_np
     u = uJ*gp(j)+(ue(1)+ue(2))/2;
     Uders=bspbasisDers(knotU,pu,u,1); 
     basis(:,j) = Uders(1,:);
     basis_grad(:,j) = Uders(2,:);
     u_jacobian(i,j) = uJ*gw(j);
    end
u_ele_basis_funcs{i} = basis;
u_ele_basis_grad{i} = basis_grad;
end

v_np = pv + t + 1;
[gp,gw]=quad_info_1d(v_np); % The quadrature points and quadrature weights in the [-1,1].
v_ele_basis_funcs = cell(vNoEs,1); 
v_ele_basis_grad  = cell(vNoEs,1); 
v_knot_span = zeros(vNoEs,1);
v_jacobian = zeros(vNoEs, v_np);
basis = zeros(pv+1,v_np);
basis_grad = basis;

for i=1:vNoEs
	ve = VBreaks(i:i+1);
    vJ=(ve(2)-ve(1))/2;
    v_knot_span(i) = findspan(knotV,pv,VBreaks(1));
    for j=1:v_np
     v = vJ*gp(j)+(ve(1)+ve(2))/2;
     Vders=bspbasisDers(knotV,pv,v,1); 
     basis(:,j) = Vders(1,:);
     basis_grad(:,j) = Vders(2,:);
     v_jacobian(i,j) = vJ*gw(j);
    end
v_ele_basis_funcs{i} = basis;
v_ele_basis_grad{i} = basis_grad;
end

w_np = pw + t + 1;
[gp,gw]=quad_info_1d(w_np); % The quadrature points and quadrature weights in the [-1,1].
w_ele_basis_funcs = cell(wNoEs,1);
w_ele_basis_grad  = cell(wNoEs,1);
w_knot_span = zeros(wNoEs,1);
w_jacobian = zeros(wNoEs, w_np);
basis = zeros(pw+1,w_np);
basis_grad = basis;

for i=1:wNoEs
	we = WBreaks(i:i+1);
    wJ=(we(2)-we(1))/2;
    w_knot_span(i) = findspan(knotW,pw,WBreaks(1));
    for j=1:w_np
     w = wJ*gp(j)+(we(1)+we(2))/2;
     Wders=bspbasisDers(knotW,pw,w,1); 
     basis(:,j) = Wders(1,:);
     basis_grad(:,j) = Wders(2,:);
     w_jacobian(i,j) = wJ*gw(j);
    end
w_ele_basis_funcs{i} = basis;
w_ele_basis_grad{i} = basis_grad;
end



n_gps = u_np*v_np*w_np; 



ele_F   =   cell(NoEs,1);    % zeros(NoEs,DIM,n_gps);
ele_DF =   cell(NoEs,n_gps);   % zeros(NoEs,n_gps,DIM,DIM);

Jacobian = zeros(NoEs,n_gps);

for k=1:wNoEs
	for j=1:vNoEs
		for i=1:uNoEs
			e = i + (j-1)*uNoEs + (k-1)*uNoEs*vNoEs;
            u_basis_i =    u_ele_basis_funcs{i};
            v_basis_j =    v_ele_basis_funcs{j};
            w_basis_k =  w_ele_basis_funcs{k};
            uv_basis_ij = kron(v_basis_j,u_basis_i);
            uvw_basis_ijk = kron(w_basis_k,uv_basis_ij);

            
            Du_basis_i =    u_ele_basis_grad{i};
            Du_v_basis_ij = kron(v_basis_j,Du_basis_i);
            Du_vw_basis_ijk = kron(w_basis_k,Du_v_basis_ij);

            
            Dv_basis_j =  v_ele_basis_grad{j};
            Dv_u_basis_ij = kron(Dv_basis_j,u_basis_i);
            Dv_uw_basis_ijk = kron(w_basis_k,Dv_u_basis_ij);

            
            
            Dw_basis_k =  w_ele_basis_grad{k};
            Dw_uv_basis_ijk = kron(Dw_basis_k,uv_basis_ij);
        
            
            uspan =  u_knot_span(i);
            vspan =  v_knot_span(j);
            wspan = w_knot_span(k);
            
            w = weights(uspan-pu:uspan,vspan-pv:vspan,wspan-pw:wspan);
            w = w(:);
            w_funcs = w'*uvw_basis_ijk; % 现在算出了权重函数在 u_np*v_np*w_np 个积分点处的�?
            
            
            P_x = ConPts(uspan-pu:uspan,vspan-pv:vspan,wspan-pw:wspan,1);
            P_x = P_x(:);
            w_P_x = (w.*P_x)';
            F_x = (w_P_x* uvw_basis_ijk)./w_funcs;
            
            
            P_y = ConPts(uspan-pu:uspan,vspan-pv:vspan,wspan-pw:wspan,2);
            P_y = P_y(:);
            w_P_y = (w.*P_y)';
            F_y = (w_P_y * uvw_basis_ijk)./w_funcs;
            
            P_z = ConPts(uspan-pu:uspan,vspan-pv:vspan,wspan-pw:wspan,3);
            P_z = P_z(:);
            w_P_z = (w.*P_z)';
            F_z = (w_P_z* uvw_basis_ijk)./w_funcs;
            
            
            ele_F{e} = [F_x;F_y;F_z];
            
            Dw_u  = w'*Du_vw_basis_ijk; 
            Dw_v  = w'*Dv_uw_basis_ijk; 
            Dw_w  = w'*Dw_uv_basis_ijk; 
            
        
            DF_x_u  = (w_P_x*Du_vw_basis_ijk - F_x.*Dw_u )./w_funcs;
            DF_x_v  = (w_P_x*Dv_uw_basis_ijk - F_x.*Dw_v )./w_funcs;
            DF_x_w = (w_P_x*Dw_uv_basis_ijk - F_x.*Dw_w )./w_funcs;

            DF_y_u  = (w_P_y*Du_vw_basis_ijk - F_y.*Dw_u )./w_funcs;
            DF_y_v  = (w_P_y*Dv_uw_basis_ijk - F_y.*Dw_v )./w_funcs;
            DF_y_w = (w_P_y*Dw_uv_basis_ijk - F_y.*Dw_w )./w_funcs;

            DF_z_u  = (w_P_z*Du_vw_basis_ijk - F_z.*Dw_u )./w_funcs;
            DF_z_v  = (w_P_z*Dv_uw_basis_ijk - F_z.*Dw_v )./w_funcs;
            DF_z_w = (w_P_z*Dw_uv_basis_ijk - F_z.*Dw_w )./w_funcs;
            
            for k1=1:n_gps
                ele_DF{e,k1} = [DF_x_u(k1)  DF_x_v(k1) DF_x_w(k1);...
                                           DF_y_u(k1)  DF_y_v(k1) DF_y_w(k1);...
                                           DF_z_u(k1)  DF_z_v(k1) DF_z_w(k1);];
            end
            
        end
    end
end




for k=1:wNoEs
	for j=1:vNoEs
		for i=1:uNoEs
			e = i + (j-1)*uNoEs + (k-1)*uNoEs*vNoEs;
            u_gw = u_jacobian(i,:);         u_gw = u_gw(:);
            v_gw = v_jacobian(j,:);         v_gw  = v_gw(:);
            w_gw = w_jacobian(k,:);      w_gw = w_gw(:);
            uv_gw = kron(v_gw,u_gw); 
            uvw_gw = kron(w_gw,uv_gw);
             for k1=1:n_gps % The loop for the quadrature points in the u,v,w-directions
                DF =ele_DF{e,k1};
                Jacobian(e,k1) = abs(det(DF))*uvw_gw(k1);
            end
        end
    end
end







 nurbsInfo.ele_F = ele_F;
 nurbsInfo.ele_DF = ele_DF;
 nurbsInfo.Jacobian = Jacobian;


end
