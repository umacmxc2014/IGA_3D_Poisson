function nurbsInfo=Iga_3d_grid(knotU,pu,knotV,pv,knotW,pw,weights,Refinement)

DIM=3;

% addpath('../NURBS/')

[Qw,Ubar,Vbar,Wbar,n_dofs]=IGAknotRefineVolume(weights,knotU,pu,knotV,pv,knotW,pw,Refinement);

nurbsInfo.Qw= Qw;
nurbsInfo.Ubar=Ubar;     % The knot vector in the u-direction.
nurbsInfo.Vbar=Vbar;     % The knot vector in the v-direction.
nurbsInfo.Wbar=Wbar;   % The knot vector in the w-direction.
nurbsInfo.n_dofs=n_dofs; % The number of DOFs in the refined grid.
nurbsInfo.pu = pu;
nurbsInfo.pv = pv;
nurbsInfo.pw = pw;

UBreaks=unique(Ubar);   % u 方向上节点向量中的断点.
VBreaks=unique(Vbar);   % v 方向上节点向量中的断点.
WBreaks=unique(Wbar); % w 方向上节点向量中的断点.

nurbsInfo.UBreaks=UBreaks;
nurbsInfo.VBreaks=VBreaks;
nurbsInfo.WBreaks=WBreaks;


N1=length(Ubar)-pu-1;  %  u 方向上基函数的个数。
N2=length(Vbar)-pv-1;  %  v 方向上基函数的个数。
N3=length(Wbar)-pw-1;%  w 方向上基函数的个数。

nurbsInfo.N1=N1;
nurbsInfo.N2=N2;
nurbsInfo.N3=N3;

uNoEs=length(UBreaks)-1;     %   u 方向上的区间数。
vNoEs=length(VBreaks)-1;      %  v 方向上的区间数。
wNoEs=length(WBreaks)-1;    %  w 方向上的区间数。
NoEs=uNoEs*vNoEs*wNoEs;   %  计算区域上的区间总数。

nurbsInfo.uNoEs=uNoEs;
nurbsInfo.vNoEs=vNoEs;
nurbsInfo.wNoEs=wNoEs;
nurbsInfo.NoEs=NoEs;

Eledof=(pu+1)*(pv+1)*(pw+1);     % 一个单元上的自由度总个数。
Element=zeros(NoEs,Eledof);         % 存储网格里每个单元上的自由度的编号.
knotSpanIndex=zeros(NoEs,DIM);  % 存储每个单元上的参数开始坐标的 knot span index. 
Coordinate=zeros(NoEs,2*DIM);     % 存储每个参数单元上的坐标起止值: $[u_i, u_{i+1}] *  [v_j, v_{j+1}]*  [w_k, w_{k+1}]$.

 for k1=1:wNoEs   % 循环w方向上的全部单元.
     for j1=1:vNoEs % 循环v方向上的全部单元.
         for i1=1:uNoEs %循环 u方向上的全部单元.
	
	e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
    % 注意，参数区域的网格单元的编号顺序是：　是先把最底部那一行单元从左到右，再从下到上来排列的.
	Coordinate(e,:)=[UBreaks(i1:i1+1),VBreaks(j1:j1+1),WBreaks(k1:k1+1)];% 存储当前单元的四个参数坐标.
    u_span=findspan(Ubar,pu,UBreaks(i1));      %当前单元e上的u方向上的节点张成区间的index,即$[u_i, u_{i+1}]$.
    v_span=findspan(Vbar,pv,VBreaks(j1));      %当前单元e上的v方向上的节点张成区间的index,即$[v_j, v_{j+1}]$.
    w_span=findspan(Wbar,pw,WBreaks(k1)); %当前单元e上的w方向上的节点张成区间的index,即$[w_k, w_{k+1}]$.
	knotSpanIndex(e,:)=[u_span,v_span,w_span];
	
    local_index = 1;
    
    for k2=(w_span-pw):w_span
        for j2=(v_span-pv):v_span
            for i2=(u_span-pu):u_span
     global_index = i2+(j2-1)*N1+(k2-1)*N1*N2;
     Element(e,local_index) = global_index;
     local_index = local_index+1;
            end
        end
    end
     % 注意，这里全局自由度的编号是先让v和w方向的index固定，把u方向上的index变化，再让w方向上的index固定，让v方向的index变化，也就是:
     % N_{1,1,1}, N_{2,1,1}, ... , N_{m,1,1}; N_{1,2,1},N_{2,2,1},..., N_{m,2,1}; ...;
     % N_{1,n,1},N_{2,n,1},..., N_{m,n,1}.
end
    end
end




nurbsInfo.Element=Element;
nurbsInfo.Coordinate=Coordinate;
nurbsInfo.knotSpanIndex=knotSpanIndex;


Element_w_0 = zeros(uNoEs*vNoEs,(pu+1)*(pv+1)); % The dofs index in the elements lying on the surface w=0;
% These elements are the 2D elements, not the elements of the 3D mesh.


for j1=1:vNoEs
    for i1=1:uNoEs
        e = i1 + (j1-1)*uNoEs; % (i,j)处的单元的全局编号为e.
        u_span=findspan(Ubar,pu,UBreaks(i1));  
        v_span=findspan(Vbar,pv,VBreaks(j1));
        local_index = 1;
   for j2=(v_span - pv):v_span
       for i2=(u_span - pu):u_span
           global_index = i2 + (j2-1)*N1;
           Element_w_0(e,local_index)=global_index;
           local_index = local_index + 1;
       end
   end
   
    end
end

nurbsInfo.Element_w_0 = Element_w_0;

bnd_ele = zeros(2*uNoEs*vNoEs+2*uNoEs*wNoEs+2*vNoEs*wNoEs,1);

index  = 1;

     k1=1;
     for j1=1:vNoEs     % 循环v方向上的全部单元.
         for i1=1:uNoEs % 循环u方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     k1=wNoEs;
     for j1=1:vNoEs     % 循环v方向上的全部单元.
         for i1=1:uNoEs % 循环u方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     
     j1=1;
     for k1=1:wNoEs   % 循环w方向上的全部单元.
         for i1=1:uNoEs % 循环u方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     j1=vNoEs;
     for k1=1:wNoEs   % 循环w方向上的全部单元.
         for i1=1:uNoEs % 循环u方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     
     i1=1;
     for k1=1:wNoEs   % 循环w方向上的全部单元.
         for j1=1:vNoEs % 循环v方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     i1=uNoEs;
     for k1=1:wNoEs   % 循环w方向上的全部单元.
         for j1=1:vNoEs % 循环v方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
   bnd_ele = unique(bnd_ele);
   
   all_ele = [1:uNoEs*vNoEs*wNoEs]';
       
   all_ele(bnd_ele) = [];
   
   interior_ele = all_ele;
   
   
   nurbsInfo.bnd_ele = bnd_ele;
   
   nurbsInfo.interior_ele = interior_ele;
     

end

