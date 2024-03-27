%% Error plot
% errorplot(x, y, va, col, alp)
% Roger Koenig

function handle = errorplot(x, y, va_max, va_min,  col, alp, linewidth)

if nargin == 5
    linewidth = 1;
end
if size(x,1) ~= 1
    x = x';
end
if size(y,1) ~= 1
    y = y';
end
if size(va_max,1) ~= 1
    va_max = va_max';
end
if size(va_min,1) ~= 1
    va_min = va_min';
end

dx = find(size(x)>1);
dy = find(size(y)>1);
dv = find(size(va_min)>1);

if isempty(x)
    x = 1:size(y,dy);
end

hold on
fill([x, flipdim(x,dx),x(1)], [y+va_max, flipdim(y,dy)-flipdim(va_min,dv),y(1)+va_max(1)], col, 'linestyle', 'none','FaceAlpha', alp)
handle = plot(x,y,'color',col,'LineWidth',linewidth);
