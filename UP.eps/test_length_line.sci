s = 'test_length_line.dat';
u=file('open',s,'unknown');
//n=300;
nst=1000;
for i=800:nst
  theta_n=i/1E6;
  t = 'lat : mylat : n_layers : ' + string(i+5);
  write(u,t);
  t = 'mol : Lin : composition : (LG)1(LA0)98(LB0)1(LA0)99(LB1)1(LA0)'+ string(i-201) +'(LB2)1';
  write(u,t);
  t = 'mol : Lin : theta : ' + string(theta_n);
  write(u,t);
  t = 'start';
  write(u,t);
  t = ' ';
  write(u,t);
end // for
file('close',u);

