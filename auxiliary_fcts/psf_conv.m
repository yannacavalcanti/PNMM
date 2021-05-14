function [X] = psf_conv(X,W,H,D,mask,kernel3D)
N = size(X,1);
GX = zeros(N,W*H*D);
GX(:,mask) = X(:,mask);
% 
% % % Imgaussfilt
% % if D==1
% %     GX = reshape(GX,N,W,H);
% %     for i=1:N
% %         GX(i,:,:) = imgaussfilt(squeeze(GX(i,:,:)),sigma);
% %     end
% %     
% %     GX = reshape(GX,N,W*H);
% % else
% %     GX = reshape(GX,N,W,H,D);
% %     
% %     for i=1:N
% %         GX(i,:,:,:) = imgaussfilt3(squeeze(GX(i,:,:,:)),sigma);
% %     end
% %     GX = reshape(GX,N,W*H*D);
% % end
% if D==1
%     GX = reshape(GX,N,W,H);
%     for i=1:N
%         %2D filt
%         GX(i,:,:) = conv2(kernel3D(1,:),kernel3D(2,:),squeeze(GX(i,:,:)),'same'); 
%         %         GX(i,:,:) = imfilter(GX(i,:,:), kernel3D, 'conv', 'circular', 'same');
%     end
%     GX = reshape(GX,N,W*H);
% else
%     GX = reshape(GX,N,W,H,D);
%     
%     for i=1:N
%         for j=1:D
%             %2D filt
%             GX(i,:,:,j) = conv2(kernel3D(1,:),kernel3D(2,:),squeeze(GX(i,:,:,j)),'same');
%         end
%         GX(i,:,:,:) = reshape(conv2(1,kernel3D(3,:),squeeze(reshape(GX(i,:,:,:),W*H,D)),'same'),W,H,D);
%              %         GX(i,:,:,:) = imfilter(GX(i,:,:,:), kernel3D, 'conv', 'circular', 'same');
%     end
%     
%     GX = reshape(GX,N,W*H*D);
% end
X(:,mask) = GX(:,mask);