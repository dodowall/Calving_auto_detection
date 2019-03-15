function [x_poly,y_poly] = mask2poly(mask,xmin,ymin)

[nrow, ncol] = size(mask);

x_poly = zeros(ncol*2,1);
y_poly = zeros(ncol*2,1);

for ic = 1:ncol
    id = find(mask(:,ic));
    if length(id>0)
        y_poly(ic) = min(id)+ymin;
        x_poly(ic) = ic+xmin;
        y_poly(end-ic+1) = max(id)+ymin;
        x_poly(end-ic+1) = ic+xmin;
    end
end;
y_poly = y_poly(x_poly>0);
x_poly = x_poly(x_poly>0);