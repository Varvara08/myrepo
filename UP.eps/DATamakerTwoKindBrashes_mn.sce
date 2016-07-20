clear;
lines(0);

np=200;
q=2;
for w=2.08

s0 = '~/imc/StarBrush/Neutral/n'+ string(np) + '/q' + string(q) + 'den/chi0/mn' + string(w) + '/' ;

nf=66;//number of the latest file
//nf=10;
i=1;
N=np*(q+1);//mass DEN VITAL!!!
for j=1:1:nf
sig=0.01*j;
//  v=j;
  s = s0 + 'CSBrush_mn' + string(w) + '_' + string(j) + '.pro';
  a=fscanfMat(s);
  [m,n]=size(a);

  zbp_av=0; //DEN
  zbpL_av=0; // Line
  zbp1L_av=0; // line mon=200
  bp_norm=0; // norm DEN
  bpL_norm=0; // norm Line
  bp1L_norm=0; // norm Line mon=200
  zbp2_av = 0; // quadrat =) DEN
  zbp2L_av = 0; // quadrat =) Line
  zbp12L_av = 0; // -/- mon =200

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
  phiL2=0;
  
  z_av=0;
  zL_av=0;
  z2_av=0;
  zL2_av=0;
  H=0;
  HL=0;
  
  for k=1:m
    b(k,1)=a(k,1); //layer
    b(k,2)=a(k,7); // phi_tot
    b(k,3)=a(k,13)/sig; // nB0 (branching point)
    b(k,4)=a(k,15)/sig/q; // nB1 (end point)
    b(k,5)=a(k,8)*1E6/(w*np); // phiL_tot
    b(k,6)=0;   //a(k,18)*1E6; // nB0L (simular branch. point; mon=100)
    b(k,7)=a(k,18)*1E6; //nB2L (end point)
//    b(k,8)=a(k,20)*1E6; //nB1 (mon=200)
                              
    zbp_av=zbp_av+k*b(k,3);
    zbp2_av=zbp2_av+k*k*b(k,3);
    bp_norm=bp_norm+b(k,3);
                              
//    zbpL_av=zbpL_av+k*b(k,6);
//    zbp2L_av=zbp2L_av+k*k*b(k,6);
//    bpL_norm=bpL_norm+b(k,6);
                                
//    zbp1L_av=zbp1L_av+k*b(k,8);
//    zbp12L_av=zbp12L_av+k*k*b(k,8);
//    bp1L_norm=bp1L_norm+b(k,8);
                                 
                              
    ze_av=ze_av+k*b(k,4);
    ze2_av=ze2_av+k*k*b(k,4);
    e_norm=e_norm+b(k,4);
                              
    zeL_av=zeL_av+k*b(k,7);
    ze2L_av=ze2L_av+k*k*b(k,7);
    eL_norm=eL_norm+b(k,7);
                              
    z_av=z_av+k*b(k,2);
    z2_av=z2_av+k*k*b(k,2);
                              
     zL_av=zL_av+k*b(k,5);
     zL2_av=zL2_av+k*k*b(k,5);
     
    phi_norm=phi_norm+b(k,2);
    phiL_norm=phiL_norm+b(k,5)
    
    phi2=phi2+b(k,2)*b(k,2);
    phi_av=phi_av+phi2/phi_norm;
    
    phiL2=phiL2+b(k,5)*b(k,5);
    phiL_av=phiL_av+phiL2/phiL_norm;
    
    end

zbp_av=zbp_av/bp_norm;
zbp2_av=zbp2_av/bp_norm;

//zbpL_av=zbpL_av/bpL_norm;
//zbp2L_av=zbp2L_av/bpL_norm;

//zbp1L_av=zbp1L_av/bp1L_norm;
//zbp12L_av=zbp12L_av/bp1L_norm;

ze_av=ze_av/e_norm;
ze2_av=ze2_av/e_norm;

zeL_av=zeL_av/eL_norm;
ze2L_av=ze2L_av/eL_norm;

z_av=z_av/phi_norm;
z2_av=z2_av/phi_norm;

zL_av=zL_av/phiL_norm;
zL2_av=zL2_av/phiL_norm;

phi_av=phi_av/phi_norm;
H=sig*N/phi_av;

res(i,1)=sig;    //sigma
//res(i,2)=layer; 
res(i,2)=zbp_av;    // точка ветвления DEN
//res(i,3)=zbp2_av;
res(i,3)=sqrt(zbp2_av-zbp_av^2);  // Дисперсия Т.В.
res(i,4)=ze_av; //  среднее положение концевой группы DEN
//res(i,6)=ze2_av;
res(i,5)=sqrt(ze2_av-ze_av^2); // дисперсия п.к.г
res(i,6)=z_av; // phi DEN
res(i,7)=sqrt(z2_av-z_av^2); //  дисперсия phi DEN
res(i,8)=zL_av; //  phi LINE
res(i,9)=sqrt(zL2_av-zL_av^2); // дисперсия phi LINE
res(i,10)=zeL_av; //  среднее положение концевой группы LINE
res(i,11)=sqrt(ze2L_av-zeL_av^2); // дисперсия п.к.г line
res(i,12)=zbpL_av; //  точка n=100 in line
res(i,13)=sqrt(zbp2L_av-zbpL_av^2); // дисперсия
res(i,14)=H; //  высота щётки
res(i,15)=zbp1L_av; //  точка n=200 in line
res(i,16)=sqrt(zbp12L_av-zbp1L_av^2); // дисперсия
res(i,17)=zeL_av/eL_norm; // норм



// resul(i,1)=sig; //sigma
// resul(i,2)=zbp_av;    // точка ветвления DEN
// resul(i,3)=z_av; // phi DEN
// resul(i,4)=zeL_av; //  среднее положение концевой группы LINE
// resul(i,5)=sqrt(ze2L_av-zeL_av^2); // дисперсия п.к.г line
// resul(i,6)=zbpL_av; //  точка n=100 in line
// resul(i,7)=sqrt(zbp2L_av-zbpL_av^2); // дисперсия

i=i+1;
  //
  sout = s0 + 'CSBrush_mn' + string(w) + '_' + string(j) + '.dat';
  u=file('open',sout,'unknown');
  write(u, b, '(f12.6,32(e16.6))');
  file('close',u);
end
  //s0 = '/home/alexk/AlexK/IMC/StarBrush/Neutral/n'+ string(np) + '/q' + string(q) + 'den/chi0.5/m' + string(w) + '/';
  sout = s0 + 'zav_q' + string(q) +'_mn'+ string(w) + '_np' + string(np) + '.dat';
  u=file('open',sout,'unknown');
  write(u, res, '(f12.6,32(e16.6))');
  file('close',u);
  
//    s0 = '/home/alexk/AlexK/IMC/StarBrush/Neutral/n=' + string(np) +  '/Result/ZLE/'+ 'q=' + string(q) + '/dat/';
//  sout = s0 + 'result_q=' + string(q) +'_m='+ string(w) + '_n=' + string(np) + '.dat';
//  u=file('open',sout,'unknown');
//  write(u, resul, '(f12.6,32(e16.6))');
//  file('close',u);
end
