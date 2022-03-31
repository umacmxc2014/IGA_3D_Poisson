function [A,rhs]=Iga_3d_bc(A,rhs,N1,N2,N3)

% The global DOF index for (i,j,k) is  i + (j-1)*N1 + (k-1)*N1*N2.

bottom_face_dof =1:N1*N2;                              %  The face  w=0，即底面。
up_face_dof        = (1:N1*N2) + (N3-1)*N1*N2; % The face  w=1，即顶面。

front_face_dof    =  zeros(1,N2*N3);                  %  The face  u=1，即正面。
local_index = 1;
for k=1:N3
    for j=1:N2
        global_index = N1 + (j-1)*N1 +  (k-1)*N1*N2;
        front_face_dof(local_index) = global_index;
        local_index = local_index + 1;
    end
end


back_face_dof    =  zeros(1,N2*N3);                  %  The face  u=0，即背面。
local_index = 1;
for k=1:N3
    for j=1:N2
        global_index = 1 + (j-1)*N1 +  (k-1)*N1*N2;
        back_face_dof(local_index) = global_index;
        local_index = local_index + 1;
    end
end

        
left_face_dof = zeros(1,N1*N3);  %  The face  v=0，即左面。
local_index = 1;
for k=1:N3
    for i=1:N1
        global_index = i +  (k-1)*N1*N2;
        left_face_dof(local_index) = global_index;
        local_index = local_index + 1;
    end
end


right_face_dof = zeros(1,N1*N3);  %  The face  v=1，即右面。
local_index = 1;
for k=1:N3
    for i=1:N1
        global_index = i +(N2-1)*N1 +  (k-1)*N1*N2;
        right_face_dof(local_index) = global_index;
        local_index = local_index + 1;
    end
end


 be=[bottom_face_dof,up_face_dof,front_face_dof,back_face_dof,left_face_dof,right_face_dof];

be=unique(be);
A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
A(be,be)=eye(length(be));
