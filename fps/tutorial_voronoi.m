% Numerical geometry of nonrigid shapes.
% (C) Alexander & Michael Bronstein, 2008.
% All rights reserved.

rand('seed', 0);

% Farthest point sampling demo
load mario shape;

% Sampling contains the indices of the points on the mesh.
N      = 50;                                  % Target sample size.
sample = round(rand*(length(shape.X)-1)) + 1; % Initialize with a random point.
D      = repmat(Inf, [length(shape.X) N]);    % Distance maps.
d      = repmat(Inf, [length(shape.X) 1]);

% Colormap
cmap = hsv(N); 
perm = randperm(N);

figure(1);
set(gcf, 'Position',[1 36 1600 1085]);
clf;
set(gcf, 'Color', [1 1 1]);
%caxis([1 N]);
colorbar;

for k = 1:N-1,

    % Compute distance map from sample on the shape.
    u = repmat(Inf, [length(shape.X) 1]);
    u(sample(end)) = 0;
    D(:,k) = fastmarch(shape.TRIV, shape.X, shape.Y, shape.Z, u, set_options('mode', 'single'));
    %d = min(D,[],2);
    d = min(d, D(:,k));
    [r, idx] = max(d);

    % Visualize Voronoi map
    [vor, edges] = voronoi_tessellation(shape, sample, D);
    trisurf(vor.TRIV, vor.X, vor.Y, vor.Z, perm(vor.tri_labels));
    hold on;
    h = line(edges.X, edges.Y, edges.Z);
    set(h, 'Color', [0 0 0]);
    h = plot3(shape.X(sample), shape.Y(sample), shape.Z(sample), 'ok');
    set(h,'MarkerSize', 5, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
    hold off;
    axis image; axis off; shading flat; lighting phong; 
    colormap(cmap);
    view([-15 25]);
    camlight head;
    caxis([1 N]);
    
    title(sprintf('N = %d   r = %-.2f', k+1, r));
    
    drawnow;
    
    sample = [sample; idx];

end




