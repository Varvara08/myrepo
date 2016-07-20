clear;
lines(0);

cutoff = 1.0;
np=50;
q=2;
chi=0;
//
nf=66;//number of the latest sigma
//
kw=0;
for w=2.1:0.1:2.9
    kw=kw+1;
    s0 = '~/imc/StarBrush/Neutral/n' + string(np) + '/q' + string(q) + 'den/chi'+string(chi)+'/mn' + string(w) + '/' ;
    dmax = 1;
    sigtr=0;
//    sig1=0;
//    sig2=0;
    nmax_prev=1;
        s = s0 + 'zav_q'+string(q)+'_mn' + string(w) + '_np' + string(np) + '.dat';
        a=fscanfMat(s);
        [m,n]=size(a);
        x1 = a(1,11); //fluct terminal group LB2 first position
        x2 = a(2,11); //same thing 2nd position
        nmax=0;
        for k=3:nf
            x3=a(k,11);
            if((x1 < x2)&(x2 > x3)&(x2 > cutoff)) then
                nmax = nmax + 1;
                xm(nmax)=k-1;
                ym(nmax)=x2;
                sigtr=a(k,1)
            end // if
            x1=x2;
            x2=x3;
            sigtr1=a(k,1)
        end // for k
        if nmax==2 then
            if(abs(ym(1)-ym(2)) < dmax) then
                dmax=abs(ym(1)-ym(2));
            end
        end
    res(kw,1)=w*np;
    res(kw,2)=sigtr;
end

s0 = '~/imc/StarBrush/Neutral/n' + string(np) + '/Result/Diagrama/spibi/';
sout = s0 + 'fluct_spibi_q' + string(q) + '_np' + string(np) + '_M_sigtr.dat';
u=file('open',sout,'unknown');
write(u, res, '(f12.6,32(e16.6))');
file('close',u);
