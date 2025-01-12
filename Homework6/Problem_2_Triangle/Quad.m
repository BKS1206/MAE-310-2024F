function val = Quad(aa, xi, eta)

% Triangle element
if aa == 1
    val = 1 - xi - eta;
elseif aa == 2
    val = xi;
elseif aa == 3
    val = eta;
else
    error('Error: value of a should be 1,2,3');
end

% Quadrilateral element
% if aa == 1
%     val = 0.25 * (1-xi) * (1-eta);
% elseif aa == 2
%     val = 0.25 * (1+xi) * (1-eta);
% elseif aa == 3
%     val = 0.25 * (1+xi) * (1+eta);
% elseif aa == 4
%     val = 0.25 * (1-xi) * (1+eta);
% else
%     error('Error: value of a should be 1,2,3, or 4.');
% end

% EOF