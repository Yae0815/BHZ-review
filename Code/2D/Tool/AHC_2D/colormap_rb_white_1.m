% @author Angus Huang rabbit4a9@gmail.com
function sky = colormap_rb_white_1
 colormap_00 = [
      0,   0, 255;
      0,   0, 255;
     10,  10, 255;
    255, 255, 255;
    255,  10,  10;
    255,   5,   5;
    255,   0,   0;
    ];


colormap_0 = colormap_00/255.0;
sky = interp1(1:size(colormap_0,1), ...
    colormap_0, ...
    linspace(1,size(colormap_0,1),1024), 'pchip');

