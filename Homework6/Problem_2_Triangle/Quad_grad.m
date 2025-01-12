function [val_xi, val_eta] = Quad_grad(aa, xi, eta)

% Triangle element
if aa == 1
    val_xi  = -1;
    val_eta = -1;
elseif aa == 2
    val_xi  =  1;
    val_eta =  0;
elseif aa == 3
    val_xi  = 0;
    val_eta = 1;
else
    error('Error: value of a should be 1,2,3');
end

% Quadrilateral element
% if aa == 1
%     val_xi  = -0.25 * (1-eta);
%     val_eta = -0.25 * (1-xi);
% elseif aa == 2
%     val_xi  =  0.25 * (1-eta);
%     val_eta = -0.25 * (1+xi);
% elseif aa == 3
%     val_xi  = 0.25 * (1+eta);
%     val_eta = 0.25 * (1+xi);
% elseif aa == 4
%     val_xi  = -0.25 * (1+eta);
%     val_eta =  0.25 * (1-xi);
% else
%     error('Error: value of a should be 1,2,3, or 4.');
% end


% EOF