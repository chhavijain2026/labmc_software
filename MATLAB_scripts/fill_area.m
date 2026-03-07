% function to draw shaded plot
% input : x : 1D array of X-values
%         Uy : 1D array of upper y-boundary
%         Ly : 1D array of lower y-boundary
%         color : color of shaded area
%         transparency : value between 1 (opaque) and 0 (transparent)
% lgentry = legend Display 

function fill_area(x,Ly,Uy,color,transparency,lgentry)

[rx,cx] = size(x);
[rly,cly] = size(Ly);
[ruy,cuy] = size(Uy);

% convert all vectors to row vector
if rx~=1
    x = x';
end;
if rly~=1
    Ly = Ly';
end;
if ruy~=1
    Uy = Uy';
end;

X = [x fliplr(x)]; % row vector of x concatenated with x vector flipped 
                   % from left to right
Y = [Uy fliplr(Ly)]; % combine Y coordinates of upper line with Y 
                     % coordinates of lower line that have been flipped
                     % from left to right to form a continuous polygon
                     
hold on
h = fill(X,Y,'b','DisplayName',lgentry);
set(h,'FaceAlpha',transparency,'EdgeColor',color,'FaceColor',color)

return;

    