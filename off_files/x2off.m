% This function uses meshlab to convert between filetypes. It's useful when
% you have to batch convert and don't want to do it manually. 
% It can be used as in this example:
% folder = 'C:\Code\nonRigidPuzzle\Data\my_scans\';
% d = dir([folder '*.obj'])
% for i=1:numel(d)
%     x2off([folder d(i).name]);
% end
% written by: Or Litnay


function x2off(filename_in)

    cmnd = ['C:' ' && ' 'cd C:\Program Files\VCG\MeshLab' ' && ' 'meshlabserver -i '...
    filename_in ' -o ' [filename_in(1:end-3) 'off'] ' -om fc vc'];
    system(cmnd);

end

