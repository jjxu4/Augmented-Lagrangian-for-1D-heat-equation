function [u,p,M,b,matrices] = approxPVEsol(params,source,endpoints,setup,a)
%Using sources find u, p, and x for setup 1=constant linear, 2=nonconstant
%linear, 3=constant nonlinear, 4=nonconstant linear find u,p, and x
%For PVE model, nodes argment is optional with default 2560. To caclulate
%adjoint, must enter nodes and a as the 6th argument.


%This section sets the default number of nodes to 2560
% if nargin >4
%   params.n = nodes;
%   params.nx=nodes;
%   params.nt=nodes;
% else
%   params.n = 2560;
% end
params.n=params.nx; %To do: change all params.n in the code.

%Change source names for adjoint
if nargin>5
    source.F=source.z1;
    source.S=source.z2;
    source.g=source.z3;
    source.psi=source.z4;
end

%Solve for x, u, p for n nodes.
x=linspace(endpoints.xstart,endpoints.xend,params.nx+1);
params.x=x;
params.h=(endpoints.xend-endpoints.xstart)/params.nx;
params.deltat=(endpoints.tend-endpoints.tstart)/params.nt;


switch setup
    case 1
        if nargin>5
            [u,p,M,b,matrices]=constlinPVE(source,params,a);
        else
            [u,p,M,b,matrices]=constlinPVE(source,params);
        end
    case 2
        if nargin>5
            [u,p,M,b,matrices]=linPVE(source,params,endpoints,a);
        else
            [u,p,M,b,matrices]=linPVE(source,params,endpoints);
        end
%     case 3
%         [u,p,M,b,matrices]=constnonlinPVE(source,params);
%     case 4
%         [u,p,M,b,matrices]=nonlinPVE(source,params,endpoints);
%         u=u';
%         p=p';
    otherwise
        disp('Enter 1 for constant linear, 2 for linear, 3 for constnat nonlinear, or 4 for nonlinear')
end

end