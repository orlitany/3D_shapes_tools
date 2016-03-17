function [surface_x1,edges] = geodesic_contour(surface_x0,ddx,th)

DX = ddx(surface_x0.TRIV);
dmin = min(DX,[],2);
dmax = max(DX,[],2);
 
t1 = find(dmin <= th & dmax >= th & DX(:,1) < DX(:,2) & DX(:,1) < DX(:,3));
t2 = find(dmin <= th & dmax >= th & DX(:,2) < DX(:,1) & DX(:,2) < DX(:,3));
t3 = find(dmin <= th & dmax >= th & DX(:,3) <= DX(:,1) & DX(:,3) <= DX(:,2));

%surface_x0.TRIV(t1,:) = surface_x0.TRIV(t1,[1 2 3]);
surface_x0.TRIV(t2,:) = surface_x0.TRIV(t2,[2 3 1]);
surface_x0.TRIV(t3,:) = surface_x0.TRIV(t3,[3 1 2]);

DX = ddx(surface_x0.TRIV);
t = find(dmin <= th & dmax >= th);
tri = surface_x0.TRIV(t,:);

d1 = DX(t,1);
d2 = DX(t,2);
d3 = DX(t,3);

lambda1 = 1-(th-d1)./(d3-d1);
lambda2 = 1-(th-d1)./(d2-d1);
lambda3 = 1-(th-d2)./(d3-d2);

v1 = [surface_x0.X(tri(:,1)).*lambda1 + surface_x0.X(tri(:,3)).*(1-lambda1) , ...
      surface_x0.Y(tri(:,1)).*lambda1 + surface_x0.Y(tri(:,3)).*(1-lambda1) , ...
      surface_x0.Z(tri(:,1)).*lambda1 + surface_x0.Z(tri(:,3)).*(1-lambda1)  ];

v2 = [surface_x0.X(tri(:,1)).*lambda2 + surface_x0.X(tri(:,2)).*(1-lambda2) , ...
      surface_x0.Y(tri(:,1)).*lambda2 + surface_x0.Y(tri(:,2)).*(1-lambda2) , ...
      surface_x0.Z(tri(:,1)).*lambda2 + surface_x0.Z(tri(:,2)).*(1-lambda2)  ];

v3 = [surface_x0.X(tri(:,2)).*lambda3 + surface_x0.X(tri(:,3)).*(1-lambda3) , ...
      surface_x0.Y(tri(:,2)).*lambda3 + surface_x0.Y(tri(:,3)).*(1-lambda3) , ...
      surface_x0.Z(tri(:,2)).*lambda3 + surface_x0.Z(tri(:,3)).*(1-lambda3)  ];

% 1<th 2,3>th
%t = find(DX(:,1) < th & DX(:,2) >= th & DX(:,3) >= th);

% d2,d3>th - v1,v2
% d3>th - v1,v3
% d2>th - v2,v3

i1 = find(d2>=th & d3>=th);
i2 = find(d1<th & d2<th);
v2(i2,:) = v3(i2,:);
i3 = find(d1<th & d3<th);
v1(i3,:) = v3(i3,:);



edges = [];
edges.X = [v1(:,1)'; v2(:,1)'];
edges.Y = [v1(:,2)'; v2(:,2)'];
edges.Z = [v1(:,3)'; v2(:,3)'];

lt = zeros(size(surface_x0.TRIV,1),1);
lt(dmax < th) = 1;
lt(dmin > th) = 0;

idx1 = length(surface_x0.X) + [1:size(v1,1)]';
idx2 = length(surface_x0.X) + size(v1,1) + [1:size(v1,1)]';
%idxt1 = size(surface_x1.TRIV,1) + [1:length(i1)];
%idxt2 = size(surface_x1.TRIV,1) + length(i1) + [1:length(i2)];
%idxt3 = size(surface_x1.TRIV,1) + length(i1) + length(i2) + [1:length(i3)];


surface_x1 = surface_x0;
surface_x1.X = [surface_x1.X(:); v1(:,1); v2(:,1)];
surface_x1.Y = [surface_x1.Y(:); v1(:,2); v2(:,2)];
surface_x1.Z = [surface_x1.Z(:); v1(:,3); v2(:,3)];

surface_x1.TRIV(t(i1),:) = [surface_x0.TRIV(t(i1),1) idx1(i1) idx2(i1)];
lt(t(i1)) = 1;

n0 = size(surface_x1.TRIV,1);
surface_x1.TRIV = [surface_x1.TRIV; [idx1(i2) idx2(i2) surface_x0.TRIV(t(i2),2)]];
surface_x1.TRIV = [surface_x1.TRIV; [idx1(i2) surface_x0.TRIV(t(i2),2) surface_x0.TRIV(t(i2),1)]];
surface_x1.TRIV = [surface_x1.TRIV; [idx1(i3) idx2(i3) surface_x0.TRIV(t(i3),3)]];
surface_x1.TRIV = [surface_x1.TRIV; [idx2(i3) surface_x0.TRIV(t(i3),1) surface_x0.TRIV(t(i3),3)]];
n1 = size(surface_x1.TRIV,1);
lt(n0+1:n1) = 1;

n0 = n1;
surface_x1.TRIV = [surface_x1.TRIV; [idx2(i1) idx1(i1) surface_x0.TRIV(t(i1),2)]];
surface_x1.TRIV = [surface_x1.TRIV; [idx1(i1) surface_x0.TRIV(t(i1),3) surface_x0.TRIV(t(i1),2)]];
surface_x1.TRIV(t(i2),:) = [surface_x0.TRIV(t(i2),3) idx2(i2) idx1(i2)];
surface_x1.TRIV(t(i3),:) = [surface_x0.TRIV(t(i3),2) idx2(i3) idx1(i3)];
n1 = size(surface_x1.TRIV,1);
lt(n0+1:n1) = 0;

lt(t(i2)) = 0;
lt(t(i3)) = 0;

surface_x1.tri_labels = lt;

l = [t(:)' [(size(surface_x0.TRIV,1)+1):size(surface_x1.TRIV,1)]];

surface_x1.TRIV(l,:) = surface_x1.TRIV(l,[3 2 1]);