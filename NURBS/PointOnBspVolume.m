function S=PointOnBspVolume(P,U,pu,u,V,pv,v,knotW,pw,w)

ndim=size(P,4);

uspan=findspan(U,pu,u);    Nu=bsplinebasis(U,pu,u);
vspan=findspan(V,pv,v);     Nv=bsplinebasis(V,pv,v);
% disp('knotW=')
% disp(knotW)
wspan=findspan(knotW,pw,w); Nw=bsplinebasis(knotW,pw,w);



S=zeros(ndim,1);

for dim=1:ndim
  temp=reshape(P(uspan-pu:uspan,vspan-pv:vspan,wspan-pw:wspan, dim),pu+1,pv+1,pw+1);
  for i=1:(pu+1)
      for j=1:(pv+1)
          for k=1:(pw+1)
              S(dim) = S(dim) + temp(i,j,k)*Nu(i)*Nv(j)*Nw(k);
          end
      end
  end
end

end
  
