
load D:\ycruzcav\Database\mat\ground_truth
I = reshape(I,128*128,64);
mask = cell(1,64);
for i=1:64
    mask{i} = find(I(:,i)>27);
end
M_img = zeros(128*128,64);
for i=1:64
M_img(mask{i},i) = 1;
end


M_img = reshape(M_img,128,128,64);
I = reshape(I,128,128,64);
for i = 1:64
   figure;
   subplot(1,2,1)
   imagesc(M_img(:,:,i));
   
   subplot(1,2,2)
   imagesc(I(:,:,i));
end


cellsize = cellfun('length',mask);
V_mask = sum(cellsize);
V_mask(1:cellsize(1)) = 128*128+mask{1};
for i=2:64
    V_mask(1+sum(cell_size(1:i-1)):cellsize(i)+sum(cell_size(1:i-1))) = 128*128*(i-1)+mask{i};
end

z = zeros(1,128*128*64);
z(V_mask) = 1; 
z = reshape(z,128,128,64);
for i = 1:64
   figure;
   imagesc(squeeze(z(:,:,i)));
end
