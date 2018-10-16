fid = fopen('zoom_2535.txt','r');
data = textscan(fid,'%d');
fclose(fid);
data_matrix = cell2mat(data);
data_2 = reshape(data_matrix,2048,2048);
image(data_2);
colormap(colorcube(256))
%colormap(jet(256))