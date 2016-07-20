clear;
lines(0);

cutoff = 1.0e-25;
np=50;
q=2;
chi=0;
//
nf=666;//number of the latest file
//
kw=0;
for w=2.06:0.02:2.08
    kw=kw+1;
    s0 = '~/imc/StarBrush/Neutral/n' + string(np) + '/q' + string(q) + 'den/chi'+string(chi)+'/mn' + string(w) + '/' ;
    dmax = 1;
    sigtr=0;
    sig1=0;
    sig2=0;
    nmax_prev=1;
    for j=1:1:nf
        sig=0.001*j;
        s = s0 + 'CSBrush_mn' + string(w) + '_' + string(j) + '.pro';
        a=fscanfMat(s);
        [m,n]=size(a);
        x1 = a(1,18); //terminal group LB2
        x2 = a(2,18); //terminal group LB2 n100 (19); n50 (18); n200 (18)
        nmax=0;
        for k=3:m
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
                sigtr=sig;
            end
        end
        if(nmax > nmax_prev) then
            sig1=sig;
        end
        if(nmax < nmax_prev) then
            sig2=sig;
        end
        nmax_prev=nmax;
    end // for j
    res(kw,1)=w*np;
    res(kw,2)=sigtr;
    res(kw,3)=sig1;
    res(kw,4)=sig2;
//end

s0 = '~/imc/StarBrush/Neutral/n' + string(np) + '/Result/Diagrama/spibi/';
sout = s0 + '2.06-2.08yam_spibi_q' + string(q) + '_np' + string(np) + '_M_sigtr_sig1_sig2.dat';
u=file('open',sout,'unknown');
write(u, res, '(f12.6,32(e16.6))');
file('close',u);
end
