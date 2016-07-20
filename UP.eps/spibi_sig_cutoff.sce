clear;
lines(0);

cutoff = 1.0e-10;
np=100;
q=2;
//
nf=101;//number of the latest file
//
kw=0;
for sig=0.1:0.1:0.5
    kw=kw+1;
    s0 = './n=' + string(np) + '/q=' + string(q) + '/sig=' + string(sig) + '/' ;
    dmax = 1;
    wtr=0;
    w1=0;
    w2=0;
    nmax_prev=1;
    for j=1:1:nf
        w=199+j;
        s = s0 + 'CSBrush_sig=' + string(sig) + '_' + string(j) + '.pro';
        a=fscanfMat(s);
        [m,n]=size(a);
        x1 = a(1,18);
        x2 = a(2,18);
        nmax=0;
        for k=3:w
            x3=a(k,18);
            if((x1 < x2)&(x2 > x3)&(x2 > cutoff)) then
                nmax = nmax + 1;
                xm(nmax)=k-1;
                ym(nmax)=x2;
            end // if
            x1=x2;
            x2=x3;
        end // for k
        if nmax==2 then
            if(abs(ym(1)-ym(2)) < dmax) then
                dmax=abs(ym(1)-ym(2));
                wtr=w;
            end
        end
        if(nmax > nmax_prev) then
            w1=w;
        end
        if(nmax < nmax_prev) then
            w2=w;
        end
        nmax_prev=nmax;
    end // for j
    res(kw,1)=sig;
    res(kw,2)=wtr;
    res(kw,3)=w1;
    res(kw,4)=w2;
end

s0 = './n=' + string(np) + '/q=' + string(q) + '/' ;
sout = s0 + 'spibi_sig_co10_q=' + string(q) + '_np=' + string(np) + '.dat';
u=file('open',sout,'unknown');
write(u, res, '(f12.6,32(e16.6))');
file('close',u);

