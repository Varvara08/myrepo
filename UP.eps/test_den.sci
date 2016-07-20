s = 'test.dat';
u=file('open',s,'unknown');
n=300;
nst=1000;
for i=1:nst
    if (i*0.6)<360 then
        theta_n=0;
    else
        theta_n=0.6*i;
    end
      t = 'mol : Den : theta : ' + string(theta_n);
      write(u,t);
      t = 'start';
      write(u,t);
      t = ' ';
      write(u,t);
end // for
file('close',u);

