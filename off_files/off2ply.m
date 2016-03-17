function off2ply(filename_in)

    cmnd = ['C:' ' && ' 'cd C:\Program Files\VCG\MeshLab' ' && ' 'meshlabserver -i '...
    filename_in ' -o ' [filename_in(1:end-3) 'ply'] ' -m fc vc'];
    system(cmnd);

end

