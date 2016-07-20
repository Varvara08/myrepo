clear;
lines(0);

np=100;
q=3;
msig = 10; // номер строки, отвечающая нужной нам плотности прививки, с которой мы будем считывать

j=0;
w=190
  s0 = '../n=' + string(np) +  '/m=' + string(w) + '/zav_q/'  ;
  s = s0 + 'zav_q=' + string(q) + '_' +'np=' + string(np) + '.dat';
  s1 = '../n=' + string(np) +  '/Result/';
  a=fscanfMat(s);
  [m,n]=size(a);

  for k=1:m
      
      for w=190:10:300
      
    j=j+1;

  b(j,1)=w; // w
  b(j,2)=a(msig,10); // zeL_av
  b(j,3)=a(k,10)
  b(j,4)=a(k,1)
  
//  b(j,k)=a(k,10)
end
end

sig = msig*0.01; // плотность прививки, соответствующая номеру строки msig
sout = s1 + 'zav_q=' + string(q) + '_np=' + string(np) + '_sig=' + string(sig) + '.dat';
u=file('open',sout,'unknown');
write(u, b, '(f12.6,32(e16.6))');
file('close',u);
