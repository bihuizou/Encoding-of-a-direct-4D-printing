function   [K,H,Pmax,Pmin,D1,D2] = surfcurvature(X,Y,Z)

% X=xx;
% Y=yy;
% Z=zz;

% First Derivatives
[Xu,Xv] = gradient(X);
[Yu,Yv] = gradient(Y);
[Zu,Zv] = gradient(Z);

% Second Derivatives
[Xuu,Xuv] = gradient(Xu);
[Yuu,Yuv] = gradient(Yu);
[Zuu,Zuv] = gradient(Zu);

[Xuv,Xvv] = gradient(Xv);
[Yuv,Yvv] = gradient(Yv);
[Zuv,Zvv] = gradient(Zv);

% Reshape 2D Arrays into Vectors
Xu = Xu(:);   Yu = Yu(:);   Zu = Zu(:); 
Xv = Xv(:);   Yv = Yv(:);   Zv = Zv(:); 
Xuu = Xuu(:); Yuu = Yuu(:); Zuu = Zuu(:); 
Xuv = Xuv(:); Yuv = Yuv(:); Zuv = Zuv(:); 
Xvv = Xvv(:); Yvv = Yvv(:); Zvv = Zvv(:); 

ru          =   [Xu Yu Zu];
rv          =   [Xv Yv Zv];
ruu         =   [Xuu Yuu Zuu];
ruv         =   [Xuv Yuv Zuv];
rvv         =   [Xvv Yvv Zvv];

% First fundamental Coeffecients of the surface (E,F,G)
E           =   dot(ru,ru,2);
F           =   dot(ru,rv,2);
G           =   dot(rv,rv,2);

m           =   cross(ru,rv,2);
p           =   sqrt(dot(m,m,2));
n           =   -m./[p p p]; 

[s,t] = size(Z);

% [nu,nv] = gradient(reshape(n,s,t,3));
% Nu = reshape(nu,[],3);
% Nv = reshape(nv,[],3);

% Second fundamental Coeffecients of the surface (L,M,N)
L           =   dot(ruu,n,2);
M           =   dot(ruv,n,2);
N           =   dot(rvv,n,2);

% Gaussian Curvature
K = (L.*N - M.^2)./(E.*G - F.^2);
K = reshape(K,s,t);

% Mean Curvature
H = (E.*N + G.*L - 2.*F.*M)./(2*(E.*G - F.^2));
H = reshape(H,s,t);

% Principal Curvatures
Pmax = H + sqrt(H.^2 - K);
Pmin = H - sqrt(H.^2 - K);

end




