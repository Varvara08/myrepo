clear;
lines(0);


np=100;
q=2;
w=240;
s0 = './n=' + string(np) + '/q=' + string(q) + '/m=' + string(w) + '/' ;

theta=1.0e-6;

nf=60;//number of the latest file
//nf=10;
i=1;
N=300;
for j=1:1:nf
  sig=0.01*j;
//  v=j;
  s = s0 + 'CSBrush_m=' + string(w) + '_' + string(j) + '.pro';
  a=fscanfMat(s);
  [m,n]=size(a);
  x1 = a(1,18);
  x2 = a(2,18);
  nmax=0;
  for k=3:m
    x3=a(k,18);
    if((x1 < x2)&(x2 > x3)) then
        nmax = nmax + 1;
        xm(nmax)=k-1;
        ym(nmax)=x2;
    end // if
    x1=x2;
    x2=x3;
  end // for k
 res(j,1)=sig;
 res(j,2)=nmax;
 for k=1:nmax
     res(j,2+2*k-1)=xm(k);
     res(j,2+2*k)=ym(k);
 end
end // for j
//  s0 = './n=' + string(np) +  '/q=' + string(q) + '/m=' + string(w) + '/zav_q/';
  sout = s0 + 'countmax_q=' + string(q) + '_np=' + string(np) + '_mp=' + string(w) + '.dat';
  u=file('open',sout,'unknown');
  write(u, res, '(f12.6,32(e16.6))');
  file('close',u);

