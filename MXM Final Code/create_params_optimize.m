function [params,endpoints] = create_params_optimize()
    %Define parameters for PDE system
    params.kref=1;
    params.phi0=0.5;
    params.lambdae=1;
    params.mue=1;
    params.muv=1;
    params.lambdav=1;
    params.delta=1;
    
    endpoints.xstart=0;
    endpoints.xend=2;
    endpoints.tstart=0;
    endpoints.tend=0.1;
    
    %Parameters needed for test case 2 source functions
    params.Uref = 1;
    params.Pref = 10;
    params.Hv = params.lambdav + 2*params.muv;
    params.Ha = params.lambdae + 2*params.mue;
    params.length=endpoints.xend-endpoints.xstart;
    params.omegax=8/(endpoints.xend-endpoints.xstart);
    params.omegat=8/(endpoints.tend-endpoints.tstart);
end