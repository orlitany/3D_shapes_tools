function shape = loadoff_color(filename)

shape = [];

f = fopen(filename, 'rt');

fgetl(f);
n = sscanf(fgetl(f), '%d %d %d');
nv = n(1);
nt = n(2);

data = fscanf(f, '%f');

shape.TRIV = reshape(data(3*nv+1:3*nv+4*nt), [4 nt])';
shape.TRIV = shape.TRIV(:,2:end) + 1;

data = data(1:7*nv);
data = reshape(data, [7 nv]);

color_per_vertex = data(4:6,:)';

shape.X = data(1,:)';
shape.Y = data(2,:)';
shape.Z = data(3,:)';

shape.color = color_per_vertex;

fclose(f);
