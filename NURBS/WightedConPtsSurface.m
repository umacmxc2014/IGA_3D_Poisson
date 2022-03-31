function Pw=WightedConPtsSurface(P,w)
%===�����NUBRS�������������ʾʱ����Ҫ�õ��Ĵ�Ȩϵ��Ŀ��Ƶ���󣻿��Ƶ����PΪm*n*ndim��
%== m��n�ֱ�ΪΪu��v�����ϵĻ���ĸ���ndimΪ���Ƶ����ڵĿռ�ά��Ϊ2����3��
%==== �����ΪPw����һ����Ȩϵ��Ŀ��Ƶ����Ϊm*n*(ndim+1)�͵ľ���


[m,n,ndim]=size(P);
Pw=zeros(m,n,ndim+1);
Pw(:,:,end)=w;
for i=1:ndim
	Pw(:,:,i)=P(:,:,i).*w;
end
end
