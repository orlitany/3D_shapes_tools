function [X, Y, Z] = noise_sphere(X, Y, Z, noise_std, min_R)
% usage example
% if (strcmp(shape.type, 'sphere'))
%     noise_std = 0.05; min_R = 0.5;
%     [shape.X, shape.Y, shape.Z] = noise_sphere(shape.X,shape.Y,shape.Z, noise_std, min_R);
%     figure; trimesh(shape.Tri, shape.X, shape.Y, shape.Z); title(['noised mesh, noise std = ' num2str(noise_std)]);
% end

[THETA,PHI,R] = cart2sph(X,Y,Z);
R = R + randn(size(R)) * noise_std;
R = max(R, min_R);
[X, Y, Z] = sph2cart(THETA,PHI,R);



