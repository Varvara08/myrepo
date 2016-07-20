clear;
lines(0);


np=200;
q=7;
sig=0.00025;
//sig1=0.01;
siglin = 1.0e-6;
s0 = '~/imc/StarBrush/Neutral/Data_17_apr/n' + string(np) + '/chi0/q' + string(q) + '/sig' + string(sig) + '/' ;

nf=601;//number of the latest file
i=1;
N=np*(q+1); // MASS VITAL
for j=1:1:nf
  s = s0 + 'CSBrush_sig' + string(sig) + '_' + string(j) + '.pro';
  a=fscanfMat(s);
  [m,n]=size(a);
  zbp_av=0;
  zbpL_av=0;
  bp_norm=0;
  bpL_norm=0;
  zbp2_av = 0;
  zbp2L_av = 0;
  
  ze_av=0;
  zeL_av=0;
  ze2_av=0;
  ze2L_av=0;
  e_norm=0;
  eL_norm=0;
  
  phi_norm=0;
  phiL_norm=0;
  phi_av=0;
  phiL_av=0;
  phi2_av=0;
  phi2L_av=0;
  phi2=0;  
  phi2L=0;
  
  z_av=0;
  zL_av=0;
  z2_av=0;
  z2L_av=0;
  H=0;
  HL=0;
  
  philtot=0.0;
  philetot=0.0;
  philbptot=0.0;
  philbp2tot=0.0;
  
  
  for k=1:m
    b(k,1)=a(k,1);
    b(k,2)=a(k,7); // phi_tot
    b(k,3)=a(k,13)/sig; // nB0 (branching point)
    b(k,4)=a(k,15)/sig/q; // nB1 (end point)
    
    b(k,5)=a(k,8); // phiL_tot
    b(k,6)=a(k,18); // nB0L (same the bp m=100)
    b(k,7)=a(k,19); // nB2L (end point)
    b(k,8)=a(k,20); // nB1L (m=200)
    
    
    zbp_av=zbp_av+k*b(k,3);
    zbp2_av=zbp2_av+k*k*b(k,3);
    bp_norm=bp_norm+b(k,3);
    
    ze_av=ze_av+k*b(k,4);
    ze2_av=ze2_av+k*k*b(k,4);
    e_norm=e_norm+b(k,4);
    
    zeL_av = zeL_av + k*b(k,7);
    ze2L_av = ze2L_av + k*k*b(k,7);
    eL_norm = eL_norm + b(k,7);
    bpL_norm = bpL_norm + a(k,18);
    
    z_av=z_av+k*b(k,2);
    z2_av=z2_av+k*k*b(k,2);
    
    phi_norm=phi_norm+b(k,2);
    
    phi2=phi2+b(k,2)*b(k,2);
    phi_av=phi_av+phi2/phi_norm;
    
    
    b(k,5)=b(k,5)/siglin; // phiL_tot
    b(k,7)=b(k,7)/siglin; // nB2L (end point)
    
    b(k,6)=b(k,6)/siglin; // nB0L (m=100) norm like end(on M)
    b(k,8)=b(k,8)/siglin; // nB1L (m=200) norm like end(on M)

    philtot=philtot+b(k,5);
    philetot=philetot+b(k,7);
    philbptot=philbptot+b(k,6);
    philbp2tot=philbptot+b(k,8);
    
    end // for k
zbp_av=zbp_av/bp_norm;
zbp2_av=zbp2_av/bp_norm;
ze_av=ze_av/e_norm;
ze2_av=ze2_av/e_norm;

zeL_av=zeL_av/eL_norm;
ze2L_av=ze2L_av/eL_norm;

z_av=z_av/phi_norm;
z2_av=z2_av/phi_norm;

phi_av=phi_av/phi_norm;
H=sig*N/phi_av;

res(i,1)=j+399; // 200 => n100, 399 => n200
res(i,2)=zbp_av;
//res(i,3)=zbp2_av;
res(i,3)=sqrt(zbp2_av-zbp_av^2);
res(i,4)=sig;
res(i,5)=ze_av;
//res(i,6)=ze2_av;
res(i,6)=sqrt(ze2_av-ze_av^2);
res(i,7)=z_av;
res(i,8)=sqrt(z2_av-z_av^2);
res(i,9)=sig;
res(i,10)=zeL_av;
res(i,11)=sqrt(ze2L_av-zeL_av^2);
res(i,12)=H;

res(i,13)=philtot;
res(i,14)=philetot;
res(i,15)=philbptot;
res(i,16)=philbp2tot;

i=i+1;
  //
  sout = s0 + 'CSBrush_sig' + string(sig) + '_' + string(j) + '.dat';
  u=file('open',sout,'unknown');
  write(u, b, '(f12.6,32(e16.6))');
  file('close',u);
end // for j
//  s0 = './n=' + string(np) +  '/q=' + string(q) + '/m=' + string(w) + '/zav_q/';
  sout = s0 + 'zav_q' + string(q) + '_np' + string(np)+ '_sig' + string(sig) + '.dat';
  u=file('open',sout,'unknown');
  write(u, res, '(f12.6,32(e16.6))');
  file('close',u);

