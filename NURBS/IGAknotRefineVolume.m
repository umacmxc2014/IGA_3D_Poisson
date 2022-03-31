function [Qw,Ubar,Vbar,Wbar,dof]=IGAknotRefineVolume(weights,knotU,pu,knotV,pv,knotW,pw,Refinement)

Ubar=knotU;Vbar=knotV; Wbar=knotW;
 
[N1,N2,N3]=size(weights); % 三个方向上的初始系数构造的三维矩阵的维数

for i=1:Refinement
	 UBreks=unique(Ubar);VBreks=unique(Vbar); WBreks=unique(Wbar);
     Xu= (UBreks(1:end-1)+UBreks(2:end))/2;
     Xv= (VBreks(1:end-1)+VBreks(2:end))/2;
     Xw=(WBreks(1:end-1)+WBreks(2:end))/2;
     temp=[Ubar,Xu]; Ubar=temp; Ubar=sort(Ubar);  
     temp=[Vbar,Xv];  Vbar=temp; Vbar=sort(Vbar); 
     temp=[Wbar,Xw];Wbar=temp;Wbar=sort(Wbar); 
end

nu=length(Ubar)-pu-1;nv=length(Vbar)-pv-1; nw=length(Wbar)-pw-1;
dof=nu*nv*nw;
Qw=zeros(nu,N2,N3);

for j=1:N2
    for k=1:N3
	[Ubar,wU]=IGAknotRefineCurve(knotU,weights(:,j,k)',pu,Refinement); % 先加密u方向上的权系数
     Qw(:,j,k)=wU';
    end
end

weights=Qw;
Qw=zeros(nu,nv,N3);

for i=1:nu
    for k=1:N3
	[Vbar,wV]=IGAknotRefineCurve(knotV,weights(i,:,k),pv,Refinement); % 再加密v方向上的权系数
     Qw(i,:,k)=wV;
end
end

weights=Qw;
Qw=zeros(nu,nv,nw);

for i=1:nu
    for j=1:nv
	[Wbar,wW]=IGAknotRefineCurve(knotW,reshape(weights(i,j,:),1,N3),pw,Refinement); % 再加密w方向上的权系数
     Qw(i,j,:)=wW;
end
end

end
