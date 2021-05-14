function plot_3D_PET(img,PET)

% figure;
hold on
x = [];
y = [];
z = [];
for i=1:PET.D
[x_tmp,y_tmp] = find(img(:,:,i)>=0.7);
z_tmp=i*ones(length(x_tmp),1);
x = [x;x_tmp];
y = [y;y_tmp];
z = [z;z_tmp];
end

scatter3(x,y,z,100,'filled','y')
% view(35,35)
axis equal;axis tight;
clear x y z

x = [];
y = [];
z = [];
for i=1:PET.D
[x_tmp,y_tmp] = find(img(:,:,i)<0.7 & img(:,:,i)>=0.4);
z_tmp=i*ones(length(x_tmp),1);
x = [x;x_tmp];
y = [y;y_tmp];
z = [z;z_tmp];
end

scatter3(x,y,z,'filled','r')
% view(35,35)
axis equal;axis tight;
clear x y z

x = [];
y = [];
z = [];
for i=1:PET.D
[x_tmp,y_tmp] = find(img(:,:,i)<0.4 & img(:,:,i)>=0.2);
z_tmp=i*ones(length(x_tmp),1);
x = [x;x_tmp];
y = [y;y_tmp];
z = [z;z_tmp];
end

scatter3(x,y,z,'filled','b')
% view(35,35)
axis equal;axis tight;
clear x y z

x = [];
y = [];
z = [];
for i=1:PET.D
[x_tmp,y_tmp] = find(img(:,:,i)<0.2);
z_tmp=i*ones(length(x_tmp),1);
x = [x;x_tmp];
y = [y;y_tmp];
z = [z;z_tmp];
end

scatter3(x,y,z,'filled','k')
view(125,45)
axis equal;axis tight;
hold off