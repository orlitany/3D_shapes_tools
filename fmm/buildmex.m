mex -output fastmarchmex my_heap.cpp unfold.cpp

load './michael10';

f = fastmarchmex('init', int32(surface.TRIV-1), double(surface.X(:)), double(surface.Y(:)), double(surface.Z(:)));

for k=1:10,
    source = repmat(Inf, [size(surface.X) 1]);
    source(round(rand*(length(surface.X)-1)+1),1) = 0;
    d = fastmarchmex('march', f, double(source));
    %d = fastmarchmex('march', f, double([k]));
    d(d>=9999999) = Inf;
    trisurf(surface.TRIV, surface.X, surface.Y, surface.Z, d(:,end)); axis image;
    drawnow;
end

fastmarchmex('deinit', f);



