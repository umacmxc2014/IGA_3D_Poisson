function [Q,wbar,knotUbar,knotVbar,knotWbar]=IGADegreeElevVolume(ConPts,weights,knotU,pu,knotV,pv,knotW,pw,t)

% Elevate the degree of B-spline basis functions to (pi +t) for i=1,2 and
% 3.


if t==0 % t=0 means that we do not elevate the degree of B-splines.
    Q=ConPts;wbar=weights;knotUbar=knotU;knotVbar=knotV; knotWbar=knotW;
end

 if t>=1
[N1,N2,N3,ndim]=size(ConPts);

UBreaks=unique(knotU);  NoUBreaks=length(UBreaks);
VBreaks=unique(knotV);  NoVBreaks=length(VBreaks);
WBreaks=unique(knotW);NoWBreaks=length(WBreaks);

     nu=N1+(NoUBreaks-1)*t;
	
     Q=zeros(nu,N2,N3,ndim);
    wbar=zeros(nu,N2,N3);
    
 for j=1:N2
     for k=1:N3
	 temp=reshape(ConPts(:,j,k,:),N1,ndim);temp=temp';
      [knotUbar,temp,wU]=DegreeElevCurve(temp,weights(:,j,k)',knotU,pu,t);temp=temp';
      Q(:,j,k,:)=temp;
      wbar(:,j,k)=wU';
 end
 end
 
ConPts=Q;weights=wbar;
nv=N2+(NoVBreaks-1)*t;
Q=zeros(nu,nv,N3,ndim);
wbar=zeros(nu,nv,N3);
 
 for i=1:nu
     for k=1:N3
	 temp=reshape(ConPts(i,:,k,:),N2,ndim);temp=temp';
      [knotVbar,temp,wV]=DegreeElevCurve(temp,weights(i,:,k),knotV,pv,t);temp=temp';
      Q(i,:,k,:)=temp;
      wbar(i,:,k)=wV;
     end
 end

 ConPts=Q;weights=wbar;
nw=N3+(NoWBreaks-1)*t;
Q=zeros(nu,nv,nw,ndim);
wbar=zeros(nu,nv,nw);

 for i=1:nu
     for j=1:nv
	 temp=reshape(ConPts(i,j,:,:),N3,ndim);temp=temp';
      [knotWbar,temp,wW]=DegreeElevCurve(temp,reshape(weights(i,j,:),1,N3),knotW,pw,t);temp=temp';
      Q(i,j,:,:)=temp;
      wbar(i,j,:)=wW;
     end
 end


    end



end

 
