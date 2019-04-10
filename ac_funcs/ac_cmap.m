function cmap = ac_cmap

map = colorblindmap/255;

cmap.target = map(2,:);
cmap.side = map(3,:);
cmap.feat = map(8,:);
cmap.dist = map(4,:);
cmap.lapse = map(1,:);