clear all;
% load your 3d point cloud and colors
load('example_scene');
scatter3(points3d(:,1),points3d(:,2),points3d(:,3),3,rgb);axis image

% To save as .off file use "saveoff_color"
saveoff_color('scene_off.off',points3d,[],rgb);

% To load an off file, use: loadoff. 
% this doesn't load the color, only xyz and TRIV (if availiable)
model_in = loadoff('model.off');

% To display meshes use showshape
showshape(model_in)


