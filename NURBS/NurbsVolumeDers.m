function [W,DW, D2W,F,DF,D2F]=NurbsVolumeDers(ConPts,knotU,knotV,knotW,weights,pu,u,pv,v,pw,w)

Uders=bspbasisDers(knotU,pu,u,2);
Nu=    Uders(1,:);
DNu=  Uders(2,:);
D2Nu=Uders(3,:);

Vders=bspbasisDers(knotV,pv,v,2);
Nv=    Vders(1,:)';
DNv=  Vders(2,:)';
D2Nv=Vders(3,:)';

Wders=bspbasisDers(knotW,pw,w,2);
Nw=    Wders(1,:)';
DNw=  Wders(2,:)';
D2Nw=Wders(3,:)';



i=findspan(knotU,pu,u);
j=findspan(knotV,pv,v);
k=findspan(knotW,pw,w);

u_index=i-pu:i;
v_index=j-pv:j;
w_index=k-pw:k;


w_ijk=weights(u_index,v_index,w_index);

P_u=ConPts(u_index,v_index,w_index,1);
P_v=ConPts(u_index,v_index,w_index,2);
P_w=ConPts(u_index,v_index,w_index,3);


DIM=3;


F=zeros(DIM,1);




W=0;
DW=zeros(1,DIM);
D2W = zeros(DIM,DIM);

for i1=1:(pu+1)
    for j1=1:(pv+1)
        for k1=1:(pw+1)
            W= W + w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*Nw(k1);
            DW(1) = DW(1) + w_ijk(i1,j1,k1)*DNu(i1)*Nv(j1)*Nw(k1);
            DW(2) = DW(2) + w_ijk(i1,j1,k1)*Nu(i1)*DNv(j1)*Nw(k1);
            DW(3) = DW(3) + w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*DNw(k1);
            
            D2W(1,1) = D2W(1,1) + w_ijk(i1,j1,k1)*D2Nu(i1)*Nv(j1)*Nw(k1);
            D2W(1,2) = D2W(1,2) + w_ijk(i1,j1,k1)*DNu(i1)*DNv(j1)*Nw(k1);
            D2W(1,3) = D2W(1,3) + w_ijk(i1,j1,k1)*DNu(i1)*Nv(j1)*DNw(k1);
            
            
            
            D2W(2,2) = D2W(2,2) + w_ijk(i1,j1,k1)*Nu(i1)*D2Nv(j1)*Nw(k1);
            D2W(2,3) = D2W(2,3) + w_ijk(i1,j1,k1)*Nu(i1)*DNv(j1)*DNw(k1);
            D2W(3,3) = D2W(3,3) + w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*D2Nw(k1);
            
            F(1)= F(1) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*Nw(k1);
            F(2)= F(2) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*Nw(k1);
            F(3)= F(3) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*Nw(k1);
            
        end
    end
end

F=F/W;


D2W(2,1) = D2W(1,2);  D2W(3,1) = D2W(1,3); D2W(3,2) = D2W(2,3);





W2 = W*W;
W3 = W*W2;


DF=zeros(DIM,DIM);

R_DNu =  (DNu*W-Nu*DW(1))/W2;    % The derivative of N_{i,pu}(u)/W(u,v,w)
R_DNv =  (DNv*W-Nv*DW(2))/W2;    % The derivative of N_{j,pv}(v)/W(u,v,w)
R_DNw = (DNw*W-Nw*DW(3))/W2; % The derivative of N_{k,pw}(w)/W(u,v,w)

D2F=zeros(DIM,DIM,DIM);
% D2F(DIM,DIM,1)存储的是 DF关于u的偏导数.
% D2F(DIM,DIM,2)存储的是 DF关于v的偏导数。
% D2F(DIM,DIM,3)存储的是 DF关于w的偏导数。


R_D2Nu =  ( (D2Nu*W-Nu*D2W(1,1))*W-(DNu*W-Nu*DW(1))*2*DW(1) )/W3; % 这里还没有把权系数 w_{ij} 放进来。
R_D2Nv =  ( (D2Nv*W-Nv*D2W(2,2))*W-(DNv*W-Nv*DW(2))*2*DW(2) )/W3;
R_D2Nw =  ( (D2Nw*W-Nw*D2W(3,3))*W-(DNw*W-Nw*DW(3))*2*DW(3) )/W3;

R_D2Nuv=  ( (DNu*DW(2) - Nu*D2W(1,2) )*W- (DNu*W - Nu*DW(1))*2*DW(2) )/W3;
R_D2Nuw=  ( (DNu*DW(3) - Nu*D2W(1,3) )*W- (DNu*W - Nu*DW(1))*2*DW(3) )/W3;
R_D2Nvw=  ( (DNv*DW(3) - Nv*D2W(2,3) )*W- (DNv*W - Nv*DW(2))*2*DW(3) )/W3;


for i1=1:(pu+1)
    for j1=1:(pv+1)
        for k1=1:(pw+1)
            DF(1,1) = DF(1,1) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNu(i1)*Nv(j1)*Nw(k1);
            DF(1,2) = DF(1,2) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNv(j1)*Nu(i1)*Nw(k1);
            DF(1,3) = DF(1,3) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNw(k1)*Nu(i1)*Nv(j1);
            
            DF(2,1) = DF(2,1) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNu(i1)*Nv(j1)*Nw(k1);
            DF(2,2) = DF(2,2) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNv(j1)*Nu(i1)*Nw(k1);
            DF(2,3) = DF(2,3) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNw(k1)*Nu(i1)*Nv(j1);
            
            DF(3,1) = DF(3,1) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNu(i1)*Nv(j1)*Nw(k1);
            DF(3,2) = DF(3,2) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNv(j1)*Nu(i1)*Nw(k1);
            DF(3,3) = DF(3,3) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*R_DNw(k1)*Nu(i1)*Nv(j1);
            
            
            D2F(1,1,1) = D2F(1,1,1) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*R_D2Nu(i1)*Nv(j1)*Nw(k1);
            D2F(1,2,1) = D2F(1,2,1) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*( R_D2Nuv(i1)*Nv(j1)+R_DNu(i1)*DNv(j1) )*Nw(k1);
            D2F(1,3,1) = D2F(1,3,1) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*( R_D2Nuw(i1)*Nw(k1)+R_DNu(i1)*DNw(k1) )*Nv(j1);
            
            D2F(2,1,1) = D2F(2,1,1) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*R_D2Nu(i1)*Nv(j1)*Nw(k1);
            D2F(2,2,1) = D2F(2,2,1) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*(R_D2Nuv(i1)*Nv(j1)+R_DNu(i1)*DNv(j1) )*Nw(k1);
            D2F(2,3,1) = D2F(2,3,1) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*(R_D2Nuw(i1)*Nw(k1)+R_DNu(i1)*DNw(k1) )*Nv(j1);
            
            D2F(3,1,1) = D2F(3,1,1) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*R_D2Nu(i1)*Nv(j1)*Nw(k1);
            D2F(3,2,1) = D2F(3,2,1) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*(R_D2Nuv(i1)*Nv(j1)+R_DNu(i1)*DNv(j1) )*Nw(k1);
            D2F(3,3,1) = D2F(3,3,1) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*(R_D2Nuw(i1)*Nw(k1)+R_DNu(i1)*DNw(k1) )*Nv(j1);
            
            
            D2F(1,2,2) = D2F(1,2,2) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*R_D2Nv(j1)*Nu(i1)*Nw(k1);
            D2F(1,3,2) = D2F(1,3,2) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*(R_D2Nvw(j1)*Nw(k1) + R_DNv(j1)*DNw(k1) )*Nu(i1);
            
            D2F(2,2,2) = D2F(2,2,2) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*R_D2Nv(j1)*Nu(i1)*Nw(k1);
            D2F(2,3,2) = D2F(2,3,2) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*(R_D2Nvw(j1)*Nw(k1) + R_DNv(j1)*DNw(k1) )*Nu(i1);
            
            D2F(3,2,2) = D2F(3,2,2) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*R_D2Nv(j1)*Nu(i1)*Nw(k1);
            D2F(3,3,2) = D2F(3,3,2) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*(R_D2Nvw(j1)*Nw(k1) + R_DNv(j1)*DNw(k1) )*Nu(i1);
            
            
             D2F(1,3,3) = D2F(1,3,3) + P_u(i1,j1,k1)*w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*R_D2Nw(k1);
             D2F(2,3,3) = D2F(2,3,3) + P_v(i1,j1,k1)*w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*R_D2Nw(k1);
             D2F(3,3,3) = D2F(3,3,3) + P_w(i1,j1,k1)*w_ijk(i1,j1,k1)*Nu(i1)*Nv(j1)*R_D2Nw(k1);


        end
    end
end


D2F(1,1,2) = D2F(1,2,1);
D2F(2,1,2) = D2F(2,2,1);
D2F(3,1,2) = D2F(3,2,1);

D2F(1,1,3) = D2F(1,3,1);
D2F(1,2,3) = D2F(1,3,2);

D2F(2,1,3) = D2F(2,3,1);
D2F(2,2,3) = D2F(2,3,2);

D2F(3,1,3) = D2F(3,3,1);
D2F(3,2,3) = D2F(3,3,2);








end
