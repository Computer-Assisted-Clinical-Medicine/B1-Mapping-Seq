function mask = block2D_mask(Nx,Ny,len_x,len_y)
    x_center = floor(Nx/2);
    y_center = floor(Ny/2);
    halflen_x = floor(len_x/2);
    halflen_y = floor(len_y/2);

    mask = zeros(Nx,Ny);
    mask(x_center-halflen_x+1:x_center+halflen_x, ...
        y_center-halflen_y+1:y_center+halflen_y) = 1;
end