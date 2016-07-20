clear;
lines(0);


np=100;
q=2;
s0 = '../n=' + string(np) +  '/qv4=' + string(q) + '/all_sig/';

nf=100;//number of the latest file
//nf=10;
i=1;

for j=1:10:nf
  sig=0.001*j;
  s = s0+'SBrush_q=' + string(q) + '_' + string(j) + '.pro';
  a=fscanfMat(s);
  [m,n]=size(a);
  zbp_av=0;
  bp_norm = 0;
  zbp2_av = 0;
  ze_av=0;
  ze2_av=0;
  e_norm=0;
  
  for k=1:m
    b(k,1)=a(k,1);
    b(k,2)=a(k,9); // phi_tot
    b(k,3)=a(k,16)/sig; // nB0 (branching point)
    b(k,4)=a(k,18)/sig/q; // nB1 (end point)
    zbp_av=zbp_av+k*b(k,3);
    zbp2_av=zbp2_av+k*k*b(k,3);
    bp_norm=bp_norm+b(k,3);
    ze_av=ze_av+k*b(k,4);
    ze2_av=ze2_av+k*k*b(k,4);
    e_norm=e_norm+b(k,4);
    end
zbp_av=zbp_av/bp_norm;
zbp2_av=zbp2_av/bp_norm;
ze_av=ze_av/e_norm;
ze2_av=ze2_av/e_norm;

res(i,1)=sig;
res(i,2)=zbp_av;
res(i,3)=zbp2_av;
res(i,4)=sqrt(zbp2_av-zbp_av^2);
res(i,5)=ze_av;
res(i,6)=ze2_av;
res(i,7)=sqrt(ze2_av-ze_av^2);
i=i+1;
  //
  sout = s0 + 'SBrush_qv4=' + string(q) + '_' + string(i) + '.dat';
  u=file('open',sout,'unknown');
  write(u, b, '(f12.6,32(e16.6))');
  file('close',u);
end
  s0 = '../n=' + string(np) +  '/qv4=' + string(q) + '/zav_q/';
  sout = s0 + 'zav_q=' + string(q) + '_np=' + string(np) + '.dat';
  u=file('open',sout,'unknown');
  write(u, res, '(f12.6,32(e16.6))');
  file('close',u);

