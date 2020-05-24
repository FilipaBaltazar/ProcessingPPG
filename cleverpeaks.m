function [yP,xP,yV,xV] = cleverpeaks(y, threshold) 
yP = [];
xP = [];
yV = [];
xV = [];

dy = diff(y);

[~, dxP] = findpeaks(-dy,'MinPeakHeight',0);

[yPaux,xPaux] = findpeaks(y);
[~,xVaux] = findpeaks(-y);
yVaux = y(xVaux);

checkP = [];
checkV = [];

for i = 1:length(dxP)
    prev_xP = find(xPaux <= dxP(i));
    if ~isempty(prev_xP)
        prev_xP = prev_xP(end);
        
        next_xV = find(dxP(i) <= xVaux);
        if ~isempty(next_xV)
            next_xV = next_xV(1);
            
            if yPaux(prev_xP) - yVaux(next_xV) >= threshold...
                    && ~ismember(prev_xP,checkP)...
                    && ~ismember(next_xV,checkV)
                yP = [yP,yPaux(prev_xP)];
                xP = [xP,xPaux(prev_xP)];
                yV = [yV,yVaux(next_xV)];
                xV = [xV,xVaux(next_xV)];
                
                
            else
                disp('false');
            end
        end
    end
    
end

    


end