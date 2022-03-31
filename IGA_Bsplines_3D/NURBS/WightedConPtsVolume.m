function Pw=WightedConPtsVolume(ConPts,weights)

[nu,nv,nw,ndim]=size(ConPts);
Pw=zeros(nu,nv,nw,ndim+1);

for i=1:nu
    for j=1:nv
        for k=1:nw
            Pw(i,j,k,ndim+1)=weights(i,j,k);
        end
    end
end



for dim=1:ndim
    
    for i=1:nu
        for j=1:nv
            for k=1:nw
	Pw(i,j,k,dim)=ConPts(i,j,k,dim)*weights(i,j,k);
            end
        end
    end
    
end


end
