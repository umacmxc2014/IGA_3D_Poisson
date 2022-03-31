function [DSu,DSv,DSw]=bspVolumeder(P,U,pu,u,V,pv,v,W,pw,w)

%　The code should be modified such that the inputs should be the values of basis functions at the Gauss quadrature points,
%    then we do not need to compute the basis functions again.
 
uspan=findspan(U,pu,u); 
vspan=findspan(V,pv,v);
wspan=findspan(W,pw,w);

ndim=size(P,4); % The dimension of physical domain is 3

temp=P(uspan-pu:uspan,vspan-pv:vspan,wspan-pw:wspan,:); 

ConPtsUbar=zeros(pu,pv+1,pw+1,ndim); 
ConPtsVbar=zeros(pu+1,pv,pw+1,ndim); 
ConPtsWbar=zeros(pu+1,pv+1,pw,ndim); 

tempU= U(uspan+1:uspan+pu)-U(uspan-pu+1:uspan);


tempV=V(vspan+1:vspan+pv)-V(vspan-pv+1:vspan);

tempW=W(wspan+1:wspan+pw)-W(wspan-pw+1:wspan);


for  dim=1:ndim
    
    for i=1:pu
        for j=1:pv+1
            for k=1:pw+1
                ConPtsUbar(i,j,k,dim) = pu*(temp(i+1,j,k,dim) -  temp(i,j,k,dim))/tempU(i);
            end
        end
    end
    
        for i=1:pu+1
           for j=1:pv
              for k=1:pw+1
                ConPtsVbar(i,j,k,dim) = pv*(temp(i,j+1,k,dim) -  temp(i,j,k,dim))/tempV(j);
            end
        end
        end
    
        for i=1:pu+1
           for j=1:pv+1
              for k=1:pw
                ConPtsWbar(i,j,k,dim) = pw*(temp(i,j,k+1,dim) -  temp(i,j,k,dim))/tempW(k);
            end
        end
        end
    

end

 
Ubar=U;
Ubar([1,end])=[];
Nu=bsplinebasis(Ubar,pu-1,u); 
Nv=bsplinebasis(V,pv,v);
Nw=bsplinebasis(W,pw,w);

DSu=zeros(ndim,1); 

for dim=1:ndim
    for i=1:pu
        for j=1:pv+1
            for k=1:pw+1
   DSu(dim)=DSu(dim) +ConPtsUbar(i,j,k,dim)*Nu(i)*Nv(j)*Nw(k);
            end
        end
    end
    
	end


Vbar=V;
Vbar([1,end])=[];
Nu=bsplinebasis(U,pu,u);
Nv=bsplinebasis(Vbar,pv-1,v);
Nw=bsplinebasis(W,pw,w);
DSv=zeros(ndim,1); 

for dim=1:ndim
 
    for i=1:pu+1
        for j=1:pv
            for k=1:pw+1
   DSv(dim)=DSv(dim) + ConPtsVbar(i,j,k,dim)*Nu(i)*Nv(j)*Nw(k);
            end
        end
    end
    
end


Wbar=W;
Wbar([1,end])=[];
Nu=bsplinebasis(U,pu,u);
Nv=bsplinebasis(V,pv,v);
Nw=bsplinebasis(Wbar,pw-1,w);

DSw=zeros(ndim,1); 

for dim=1:ndim
   
    for i=1:pu+1
        for j=1:pv+1
            for k=1:pw
   DSw(dim)=DSw(dim) + ConPtsWbar(i,j,k,dim)*Nu(i)*Nv(j)*Nw(k);
            end
        end
    end
    
end


    

end


