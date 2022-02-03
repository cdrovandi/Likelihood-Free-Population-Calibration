function fill_between(x,y1,y2,varargin)
        
    if nargin == 4
        alpha = varargin{1};
    else
        alpha = 0.2;
    end
   
    X = [x,fliplr(x)];
    Y = [y1,fliplr(y2)];
    patch(X,Y,'k','FaceAlpha',alpha,'EdgeColor','none');

end