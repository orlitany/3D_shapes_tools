function plot_mesh(N,h)

    if nargin == 1, h = zeros(size(N.VERT,1),1); end
    trisurf(N.TRIV,N.VERT(:,1),N.VERT(:,2),N.VERT(:,3),h),
    axis equal
	xlabel('X')
	ylabel('Y')
	zlabel('Z')
    rotate3d on
end
