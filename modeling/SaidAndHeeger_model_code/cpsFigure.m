function h = cpsFigure(width,height)
%cpsFigure(widthscale, heightscale)

h = figure;
Position = get(h,'Position');
Position(3) = width*Position(3);
Position(4) = height*Position(4);
set(h,'Position', Position)