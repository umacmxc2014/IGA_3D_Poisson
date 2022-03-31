function S=PointOnBspSurface(P,U,p,u,V,q,v)
%%==============����B��������������棨x(u,v),y(u,v),z(u,v)��= S(u,v)�ڲ���㣨u��v�����ĺ���ֵ��
%--- pΪu�����B��������Ĵ���qΪv�����B��������Ĵ���
% m=length(U)-p-1;%==== x����Ľڵ�����U�Ļ������Ϊm;
% n=length(V)-q-1;%====== y����Ľڵ������Ļ������Ϊn ��
%========= ����PΪ���Ƶ㹹�ɵľ����� m * n *ndim�͵Ŀ����� 
%%
ndim=size(P,3);
uspan=findspan(U,p,u);Nu=bsplinebasis(U,p,u);
vspan=findspan(V,q,v);Nv=bsplinebasis(V,q,v);

S=zeros(ndim,1);

for i=1:ndim
  temp=reshape(P(uspan-p:uspan,vspan-q:vspan,i),p+1,q+1);
 S(i)=Nu'*temp*Nv;
end
%==�����?���ɼ�PointOnBspSurface1��PointOnBspSurface2��������һ�µģ�ֻ���ڴ�������ϣ������Ե���ȻһЩ��
% S=zeros(1,1,ndim);
% for i=0:p
% for j=0:q
%    row=uspan-p+i;column=vspan-q+j;
%    S=S + P(row,column,:)*Nu(i+1)*Nv(j+1);
% end
% end
% S=reshape(S,ndim,1);



%%
%%==========================Test==================
% ���� The NuRBS Book Page 116
% U=[0 0 0 1/2 1 1 1];p=2;u=3/10;m=4;
% V=[0 0 0 1 1 1];q=2;v=6/10;n=3;
% P=zeros(m,n,3);
% P(:,:,1)=[0 3 6 9;0 3 6 9;0 3 6 9]';
% P(:,:,2)=[0 0 0 0;2 2 2 2;4 4 4 4]';
% P(:,:,3)=[0 3 3 0;2 5 5 2;0 3 3 0]';
% S=PointOnBspSurface(P,U,p,u,V,q,v)
% ========= ������Ϊ��
% S=
 %  3.060000000000000
 %  2.400000000000000
 %  3.480000000000000
 %%===========================================
 %% function S=PointOnBspSurface1(P,U,p,u,V,q,v)

% ��P=reshape(P,m*n,ndim)���������������㷨ʵ�֣�
%% �����ǵ��ѿ��Ƶ���루m*n��*ndim�͵ľ����У�����Ϊ��ԪB��������Ļ���ĵĸ���Ϊm*n����ndimΪ�ռ��ά��
% m=length(U)-p-1;%------ �ɽڵ�����U�й��������B��������ĸ���
% n=length(V)-q-1;%----- �ɽڵ�����V�������B��������ĸ������Ԫ������ܸ���Ϊm*n����
%% ndim=size(P,2);%  ------- ���Ƶ����ڿռ��ά��
%% uspan=findspan(U,p,u);Nu=bsplinebasis(U,p,u); 
%% vspan=findspan(V,q,v);Nv=bsplinebasis(V,q,v);
%% S=zeros(1,ndim);
%% for j=vspan-q:vspan
%% index=(j-1)*m;  
%% for i=uspan-p:uspan
%% S=S+P(index+i,:)*Nu(i-uspan+p+1)*Nv(j-vspan+q+1);
%% end
%% end
%% S=S';
%% end
 
 %% Test=====================
 % The Nurbs Book P116 3.8.
 % P=[0,0,0;3,0,3;6,0,3;9,0,0;0,2,2;3,2,5;6,2,5;9,2,2;0,4,0;3,4,3;6,4,3;9,4,0];
% U=[0 0 0 1/2 1 1 1];V=[0 0 0  1 1 1];
% p=2;q=2;
% u=3/10;v=6/10;
% S=PointOnBspSurface1(P,U,p,u,V,q,v)
% S =
%   3.060000000000000
%   2.400000000000000
%   3.480000000000000
 %%==========================
 
 %% ��ݶԾ����������Ĳ�ͬ������������һ��С����
%% function S=PointOnBspSurface2(P,U,p,u,V,q,v)
%%==============
%--- pΪu�����B��������Ĵ���qΪv�����B��������Ĵ���
% m=length(U)-p-1;%==== x����Ľڵ�����U�Ļ������Ϊm;
% n=length(V)-q-1;%====== y����Ľڵ������Ļ������Ϊn ��
%========= ����PΪ���Ƶ㹹�ɵľ����� m * n *ndim�͵Ŀ����� 

%% ndim=size(P,3);
%% uspan=findspan(U,p,u);Nu=bsplinebasis(U,p,u);
%% vspan=findspan(V,q,v);Nv=bsplinebasis(V,q,v);
%% S=zeros(ndim,1);
%%  for i=1:ndim
%%  temp=reshape(P(uspan-p:uspan,vspan-q:vspan,i),p+1,q+1);
%% S(i)=Nu'*temp*Nv;
%%  end
%%  end

