s = 'test.dat';
u=file('open',s,'unknown');
n=800;
s=2;
f=3;
nst=1000;
for i=1:nst
  theta_n=1*i;
  s=s+1;
  f=f+1;
  t = 'msig' + string(s) +' ='+string (f) ;
//  write(u,t);
//  t = 'start';
//  write(u,t);
//  t = ' ';
  write(u,t);
end // for
file('close',u);

//b(4,j)=a(msig2,10);
