clear;
lines(0);

function find_H(w,sig)
global np q chi cutoff hz
/////////////////////////////////////////////
s0 = '~/imc/StarBrush/Neutral/n' + string(np) + '/q' + string(q) + 'den/chi'+string(chi)+'/mn' + string(w) + '/' ;
    dmax = 1;
    nmax_prev=1;
        s = s0 + 'CSBrush_mn' + string(w) + '_' + string(sig) + '.dat';
        a=fscanfMat(s);
        [m,n]=size(a);
        x1 = a(1,2); //fluct terminal group LB2 first position
        x2 = a(2,2); //same thing 2nd position
        nmax=0;
        for k=3:m
            x3=a(k,2);
            if((x1 > x2)&(x2 > x3)&(x2 > cutoff)) then
                nmax = nmax + 1;
                hz=k
            end // if
            x1=x2;
            x2=x3;
        end // for 
////////////////////////////
endfunction

global hz;
cutoff = 1.0e-5;
np=200;
q=7;
chi=0;
//
s0 = '~/imc/StarBrush/Neutral/n'+ string(np)+'/Result/Diagrama/spibi/'
//s0 = '~/imc/GNUPLOT/compn/diag/q2/'
s = s0+'test_yam_spibi_q' + string(q) + '_np' + string(np) + '_M_sigtr_sig1_sig2.dat';
a=fscanfMat(s);
[m,n]=size(a);
kw=0;
for l=1:m
    kw=kw+1;
    w = a(l,1)/np
    sig = a(l,2)*100
    find_H(w,sig)
//end
    res(kw,1)=w*np;
    res(kw,2)=hz;
end // end for for

s0 = '~/imc/StarBrush/Neutral/n' + string(np) + '/Result/Diagrama/spibi/';
sout = s0 + 'M_H_q' + string(q) + '_np' + string(np) + '.dat';
u=file('open',sout,'unknown');
write(u, res, '(f12.6,32(e16.6))');
file('close',u);
