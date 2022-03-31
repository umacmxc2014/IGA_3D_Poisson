function [S,DF,W,DW]=NurbsVolume(ConPts,weights,knotU,pu,u,knotV,pv,v,knotW,pw,w)


ndim=size(ConPts,4);

Pw=WightedConPtsVolume(ConPts,weights);
Sw=PointOnBspVolume(Pw,knotU,pu,u,knotV,pv,v,knotW,pw,w);


S=Project(Sw);
W=Sw(end); 

[DSwu,DSwv,DSww]=bspVolumeder(Pw,knotU,pu,u,knotV,pv,v,knotW,pw,w); 
DAu=DSwu(1:ndim);DAv=DSwv(1:ndim);DAw=DSww(1:ndim);
DWu=DSwu(end);   DWv=DSwv(end);    DWw=DSww(end);
DSu=(DAu-DWu*S)/W; 
DSv=(DAv-DWv*S)/W; 
DSw=(DAw-DWw*S)/W; 
DF=[DSu,DSv,DSw];
DW=[DWu,DWv,DWw];

end
 