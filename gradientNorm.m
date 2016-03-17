function [A,Dx,Dy,At,TNEIGH] = gradientNorm(shape,rho)

x = shape.X(shape.TRIV);
y = shape.Y(shape.TRIV);
z = shape.Z(shape.TRIV);
v1 = [x(:,1)-x(:,3) y(:,1)-y(:,3) z(:,1)-z(:,3)];
v2 = [x(:,2)-x(:,3) y(:,2)-y(:,3) z(:,2)-z(:,3)];

Nt = cross(v1,v2);
At = sqrt(sum(Nt.^2,2));

a=((v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2));
b=((v1(:,1).*v2(:,1) + v1(:,2).*v2(:,2) + v1(:,3).*v2(:,3)));
c=((v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2));

d=(a.*c - b.^2).^-1;
a_= c.*d;
b_ = -b.*d;
c_ = a.*d;

s=sqrt(a_.*c_ - b_.^2);
tau=a_+c_;
t=sqrt(tau+2*s);
R_a = (t.^-1).*(a_+s);
R_b = (t.^-1).*(b_);
R_c = (t.^-1).*(c_+s);

Dodd=sparse([[1:size(shape.TRIV,1)]';[1:size(shape.TRIV,1)]'],[shape.TRIV(:,3);shape.TRIV(:,1)],[-1*ones(size(shape.TRIV,1),1) ; ones(size(shape.TRIV,1),1)],size(shape.TRIV,1),size(shape.X,1));

Deven = sparse([[1:size(shape.TRIV,1)]';[1:size(shape.TRIV,1)]'],[shape.TRIV(:,3);shape.TRIV(:,2)],[-1*ones(size(shape.TRIV,1),1) ; ones(size(shape.TRIV,1),1)],size(shape.TRIV,1),size(shape.X,1));

Dx=spdiags(R_a,0,size(R_a,1),size(R_a,1))*Dodd...
    +spdiags(R_b,0,size(R_b,1),size(R_b,1))*Deven;
Dy=spdiags(R_b,0,size(R_b,1),size(R_b,1))*Dodd...
    +spdiags(R_c,0,size(R_c,1),size(R_c,1))*Deven;
    
%calc rho per triangle - average over vertexes
id1 = [shape.TRIV(:,1), [1:size(shape.TRIV,1)]'];
id2 = [shape.TRIV(:,2), [1:size(shape.TRIV,1)]'];
id3 = [shape.TRIV(:,3), [1:size(shape.TRIV,1)]'];
TNEIGH = sparse([id1(:,1);id2(:,1);id3(:,1)],[id1(:,2);id2(:,2);id3(:,2)],ones(3*size(id1,1),1));

if ~isempty(rho)
	TriRho = (TNEIGH')*rho/3;
	betta = At.*TriRho.^2;
	betta = spdiags(betta,0,size(betta,1),size(betta,1));
	A = Dx'*betta*Dx + Dy'*betta*Dy;
else
	TriRho = (TNEIGH')*ones(size(shape.X,1),1)/3;
	betta = sqrt(At).*TriRho;
	betta = spdiags(betta,0,size(betta,1),size(betta,1));
	A = [betta*Dx;betta*Dy];

end

end