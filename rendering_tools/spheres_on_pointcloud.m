function spheres_on_pointcloud(filename,sphere_rad,refine_level)

% When you need to render pointclouds in programs like Blender
% 
% spheres_on_pointcloud(filename,sphere_rad,refine_level)
% 
% it's often useful to replace each point with a sphere. This script
% % does exactly that. 
% 
% Inputs:
% filename: Enter the pointcloud filname.off% 
% sphere_rad: Choose a radius to the spheres
% refine_level: Choose how fine you want your sphere to be.
% Notice: the finer the spheres, the heavier the file gets. 
% 
% written by Or Litany (orlitany <at> gmail <dot> com )

addpath(genpath('./../off files/'));
S = loadoff_color(filename)
shape_sphere = sphere_tri('ico',refine_level,sphere_rad);

offset = 0;
for i=1:size(S.X)
    shape_sphere_ = shape_sphere;
    shape_sphere_.X = shape_sphere.X + S.X(i);
    shape_sphere_.Y = shape_sphere.Y + S.Y(i);
    shape_sphere_.Z = shape_sphere.Z + S.Z(i);
    shape_sphere_.TRIV = shape_sphere.TRIV + offset;
    shape_sphere_.color = repmat(S.color(i,:),size(shape_sphere_.X,1),1);
    
    offset = offset+size(shape_sphere_.X,1);
    spheres(i) = shape_sphere_;
end

new_mesh.X = vertcat(spheres.X);
new_mesh.Y = vertcat(spheres.Y);
new_mesh.Z = vertcat(spheres.Z);
new_mesh.TRIV = vertcat(spheres.TRIV);
new_mesh.color = vertcat(spheres.color);
filename_out = [filename(1:end-4) '_out.off'];
saveoff_color(filename_out,[new_mesh.X new_mesh.Y new_mesh.Z],new_mesh.TRIV,new_mesh.color)

end