function [F,S,g,psi,exactu,exactp,exactut]=ARMA461_new(params,xend)
%Functions representations of F, S, g, psi, exactu, exactp, and exactut
    ta=-0.25;
    tb=0.25;
    xa=-0.25;
    xb=0.25;
    q=5;
    Fbar=0.1;
    G=@(t)heaviside(t - ta-1)-heaviside(t - tb-1);
    H=@(y,t)(heaviside(y - xa)-heaviside(y-xb)).*G(t);
    Q=@(y,t)2*(tanh(q*(t-ta-1))- tanh(q*(t-tb-1)));

    Amplitude = 1;
    %F=@(y,t)Fbar*H(y,t);
    F=@(y,t)Fbar*Q(y,t);
    S=@(y,t)0;
    g=@(t)0;
    psi=@(t)0;
    exactu=@(y,t) Amplitude * H(y,t);
    exactp=@(y,t) Amplitude * H(y,t);
    exactut=@(y,t)H(y,t);
end