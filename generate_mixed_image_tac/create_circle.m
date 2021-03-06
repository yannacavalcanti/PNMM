function circlePixels = create_circle(W,H,radius)
% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = W;
imageSizeY = H;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = 25;
centerY = 25;
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
% circlePixels is a 2D "logical" array.
% Now, display it.
% image(circlePixels) ;
% colormap([0 0 0; 1 1 1]);
% title('Binary image of a circle');