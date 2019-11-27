load dryer2

%vv = wgn(1,10,29);
%var(vv)
%exit;
%return;
N = 1000;
%N = 5000;
msd = zeros(N,1);
emse = zeros(N,1);
nn = zeros(N,1);
%x = sqrt(1/2) * (randn(N,1) + 1j*randn(N,1));
%x = sqrt(1)*(rand(N,1));

vv = sqrt(10000/24) * (randn(N,1) + 1j*randn(N,1));
v1 = sqrt(1/2) * (randn(N,1) + 1j*randn(N,1));

v2 = -1 + (2)*rand(N,1);

rayleigh = 1;
expd = 2;
v3 = abs(rayleigh*randn(N,1) + 1i*rayleigh*randn(N,1)); 
pd = makedist('Exponential','mu',0.2);

%v4 = random(pd,N,1); 
v4 = exprnd(2,[N,1]);
Pr = 0.01;
bernp1 = rand(N,1) <= Pr;
v1 = real(v1);
v2 = real(v2);
v3 = real(v3);
v4 = real(v4);
vv = real(vv);

A = vv.*bernp1;
A = real(A);
L = 10;
v = A + v4;
v = real(v);
%h = randsrc(16,1,[])


r = -1 + (1+1)*rand(1,9);

%disp(r);
S = 10;
h = randsrc(L,1,[r,0;(S*0.01) (S*0.01) (S*0.01) (S*0.01) (S*0.01) (S*0.01) (S*0.01) (S*0.01) (S*0.01) (1-(S*0.01*9))]);
%disp(h);

g = zeros(1,6000);
e = zeros(1,6000);
d = zeros(1,6000);
y = zeros(1,6000);
sieg = zeros(1,6000);

W = zeros(L,1);

Nw = L;

mew = 0.8;

ronnc = 0.001;

epsilonnnc = 20;

zeta = 0.001;

sigma = 0.05;
cook = 1;
%olpa = [1 2 3 4 5 6];
%ol = olpa(1:4);
%ol = flip(ol);
%disp(ol);
lamnnc = ronnc/mew;
%imi = 
gamnnc = 0;
p = ones(L,1);
for n=L:N
    xc = x((n-(L-1)):n);
    flip(xc);
    y(n) = (W')*(xc);
    d(n) = (h')*xc + v(n);
    
%disp("d");
    %disp(d(n));
    e(n) = d(n) - y(n);
    %disp("error")
    %disp(e)
    %disp("error")
    O = sort(flip(abs(e((n-Nw+1):n))));
    K = (2*(1+real(ceil(L*Pr))));
    %disp("K");
    %disp(K);
    %disp("Nw");
    %disp(Nw);
    %Tw = diag(ones(1,Nw-K),zeros(1,K));
    AA = ones(1,Nw-K);
    BB = zeros(1,K);
    CC = [AA,BB];
    Tw = diag(CC);
    sige = sqrt(((O)*Tw*(O'))/(Nw - K));
    %disp(sige);
    gamma = 0.5;
    %disp("error")
    %disp(e(n));
    %mew = 0 ;
    if (abs(e(n))>gamma)
        mew = (1 - (gamma/abs(e(n))));
        
    end
    
    eup = ((sqrt(2*3.14)*sige)/mew);
    rou = (5/L);
    g1 = 0.001;
    alpha = zeros(1,L);
    q = zeros(1,L);
    for cc=1:L
        
    %alpha(cc) = 
    kq = max([(rou*max(g1,abs(W'))),abs(W(cc))]);  
    %disp(W(cc));
    %disp(kq);
    alpha(cc) = kq;
    %break;
    end
    ss = sum(alpha);
    for cc=1:L
        %break;
        q(cc) = (alpha(cc)/ss);
        
    end
    Q = diag(q);
    mx = min((e(n)*e(n)),eup);
    %g = mean(abs(W));
    kj = ((W.*p)/(sum(p)));
    g = kj;
    f = (sign(g-abs(W)) + 1)/2;
    F = diag(f);
    fe = sign(abs(e(n)))/(((xc')*Q*(xc))+zeta);
    %fe = sign(abs(e(n)))/(((xc')*(xc))+zeta);
    g = ((-1*F*sign(W))./(1 + (epsilonnnc*abs(W))));
    %disp(size(f));
    %disp(size(Q));
    %disp(size(xc));
    W = ((W +  (mew*fe*Q*(xc))) + (mew*lamnnc*g));
    %W = (W + (mew*fe*(xc)) + (mew*lamnnc*g));
    %break;
    MSD = -1*10*log10(norm(h-W)*norm(h-W));
    %EMSE = 10*log10((((h-W)')*(xc))*(((h-W)'*(xc))));
    msd(n) = MSD;
    %emse(n) = EMSE;
    %if (n%10 == 0)
    nn(n) = n;
    %end
    
    kj = ((W.*p)/sum(p));
    %disp(size(kj));

    p = ((sign(W-(kj))+1)./2);
    if (cook<10)
     %   disp(p);
    end
    
end
%plot(nn,msd,nn,emse);
plot(nn(L:7000),(msd(L:7000)),'*');
%MSD = 10*log10((h-W));
%EMSE = 
