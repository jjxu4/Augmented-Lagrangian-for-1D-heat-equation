clear all
%Define parameters for PDE system
params.kref=1;
params.phi0=0.5;
params.lambdae=1;
params.mue=1;
params.muv=0.25;
params.lambdav=0.0774;
params.delta=1;

endpoints.xstart=0;
endpoints.xend=1;
endpoints.tstart=0;
endpoints.tend=0.1;


%Parameters needed for test case 2 source functions
params.Uref = 0.1;
params.Pref = 1;
params.Hv = params.lambdav + 2*params.muv;
params.Ha = params.lambdae + 2*params.mue;
params.length=endpoints.xend-endpoints.xstart;
params.omegax=8/(endpoints.xend-endpoints.xstart);
params.omegat=8/(endpoints.tend-endpoints.tstart);


y=[5,10,20,40,80,160,320];
%y=[320];
erroru=zeros(length(y),1);
errorp=zeros(length(y),1);
for k=1:length(y)
    n=y(k); %2560 is the default if it is not specified for approxPVEsol
    params.nx=n;
    params.nt=n;
    params.deltat=(endpoints.tend-endpoints.tstart)/n;

    [F,source.S,source.g,source.psi,exactu,exactp]=tc2source(params,endpoints.xend);
    source.F=zeros(n+1,n);
    eta=(endpoints.xend-endpoints.xstart)/n;
    
    %source.S=@(x,t)(x-x).*(t-t);
    for i=1:n+1
        ix=endpoints.xstart+(i-1)*eta;
        for j=1:n
            it=endpoints.tstart+j*params.deltat;
            source.F(i,j)=F(ix,it);
        end
    end
    %boundaryconditions for xstart
    %Should be a constant for constant cases and a function for time dependent
    %cases
    source.bcu=@(y)y-y;
    source.bcp=@(y)y-y;

    %Initial Condition
    source.IC=@(y)y-y;

    %Set parameter to describe case. 1=constantlinear, 2=linear,
    %3=constantnonlinear, 4=nonlinear.
    setup=2;



    %x,u,p approximated. Can add 6th entry to change 2560 to a different numbers
    [u,p,M,b,matrices]=approxPVEsol(params,source,endpoints,setup,n);


    [erroru(k),errorp(k)]=L2Error(exactu,exactp, u',p',endpoints);
end
h=(endpoints.xend-endpoints.xstart)./y;
figure 
fig1=tiledlayout(1,2);
nexttile
loglog(h,erroru,'o')
xlabel('h')
ylabel('$\|u-u_h\|_{L^2(0,T;L^2(0,1))}$', 'interpreter', 'latex')

nexttile
loglog(h,errorp,'o')
xlabel('h')
ylabel('$\|p-p_h\|_{L^2(0,T;L^2(0,1))}$','interpreter','latex')

