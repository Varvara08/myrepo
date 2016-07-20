clear;
lines(0);


np=100;
q=2;

for w=190

s0 = '../n=' + string(np) +  '/q=' + string(q) + 'den'  +  '/m=' + string(w) + '/' ;

nf=670;//number of the latest file
//nf=10;
i=1;
N=300;
for j=1:1:nf
sig=0.001*j;
//  v=j;
  s = s0 + 'CSBrush_m=' + string(w) + '_' + string(j) + '.pro';
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
  
  for k=1:m
    b(k,1)=a(k,1); //layer
    b(k,2)=a(k,7); // phi_tot
    b(k,3)=a(k,13)/sig; // nB0 (branching point)
    b(k,4)=a(k,15)/sig/q; // nB1 (end point)
    
    b(k,5)=a(k,8)*1E5/w; // phiL_tot
    b(k,6)=a(k,18)*1E5; // nB0L (simular branch. point; mon=100)
    b(k,7)=a(k,19)*1E5; //nB2L (end point)
                              
    zbp_av=zbp_av+k*b(k,3);
    zbp2_av=zbp2_av+k*k*b(k,3);
    bp_norm=bp_norm+b(k,3);
                              
    zbpL_av=zbpL_av+k*b(k,6);
    zbp2L_av=zbp2L_av+k*k*b(k,6);
    bpL_norm=bpL_norm+b(k,6);
                              
    ze_av=ze_av+k*b(k,4);
    ze2_av=ze2_av+k*k*b(k,4);
    e_norm=e_norm+b(k,4);
                              
    zeL_av=zeL_av+k*b(k,7);
    ze2L_av=ze2L_av+k*k*b(k,7);
    eL_norm=eL_norm+b(k,7);
                              
    z_av=z_av+k*b(k,2);
    z2_av=z2_av+k*k*b(k,2);
                              
    phi_norm=phi_norm+b(k,2);
    
    phi2=phi2+b(k,2)*b(k,2);
    phi_av=phi_av+phi2/phi_norm;
    
    
    
    end
zbp_av=zbp_av/bp_norm;
zbp2_av=zbp2_av/bp_norm;

zbpL_av=zbpL_av/bpL_norm;
zbp2L_av=zbp2L_av/bpL_norm;

ze_av=ze_av/e_norm;
ze2_av=ze2_av/e_norm;

zeL_av=zeL_av/eL_norm;
ze2L_av=ze2L_av/eL_norm;

z_av=z_av/phi_norm;
z2_av=z2_av/phi_norm;

phi_av=phi_av/phi_norm;
H=sig*N/phi_av;

res(i,1)=sig;
res(i,2)=zbp_av; //DEN точка ветвления
//res(i,3)=zbp2_av;
res(i,3)=sqrt(zbp2_av-zbp_av^2);
res(i,4)=ze_av; //DEN концевая группа
//res(i,6)=ze2_av;
res(i,5)=sqrt(ze2_av-ze_av^2);
res(i,6)=z_av; 
res(i,7)=sqrt(z2_av-z_av^2);
res(i,8)=zeL_av; //Line концевая группа
res(i,9)=sqrt(ze2L_av-zeL_av^2);
res(i,10)=zbpL_av; //Line mon=100
res(i,11)=sqrt(zbp2L_av-zbpL_av^2);
res(i,12)=H;
res(i,13)=zeL_av/eL_norm;

i=i+1;
  //
  sout = s0 + 'CSBrush_m=' + string(w) + '_' + string(j) + '.dat';
  u=file('open',sout,'unknown');
  write(u, b, '(f12.6,32(e16.6))');
  file('close',u);
end
  s0 = '../n=' + string(np) +  '/q=' + string(q) + 'den'   + '/m=' + string(w) + '/zav_q/';
  sout = s0 + 'zav_q=' + string(q) +'_m='+ string(w) + '_np1111=' + string(np) + '.dat';
  u=file('open',sout,'unknown');
  write(u, res, '(f12.6,32(e16.6))');
  file('close',u);
end
