clear;
clear gt gf;
lines(0);

// --------------------------------------------------------------------------

function f = calc_propagators(w,d)
global gt gf np lambda;
//
// Array initialization with zeros
// gt
for s=1:np
  for z=1:np+1
    gt(s,z)=0.0;
  end // for z
end // for s
// gf
for s=1:np
  for z=1:np+d+1
    gf(s,z)=1.0;
  end // for z
end // for s
//
// Initial condition
gt(1,1)=w;
for z=1:d
  gf(1,z)=w;
end // for z

// Start recurrence
for s=2:np
  gt(s,1)=(4.0*gt(s-1,1)+gt(s-1,2))*w;
  gf(s,1)=(4.0*gf(s-1,1)+gf(s-1,2))*w;
  for z=2:s
    if z<=d then
      wx=w;
    else
      wx=1.0;
    end // if    
    gt(s,z)=(gt(s-1,z-1)+4.0*gt(s-1,z)+gt(s-1,z+1))*wx;
    gf(s,z)=(gf(s-1,z-1)+4.0*gf(s-1,z)+gf(s-1,z+1))*wx;
  end // for z
  for z=s+1:np+d
    if z<=d then
      wx=w;
    else
      wx=1.0;
    end // if    
    gf(s,z)=(gf(s-1,z-1)+4.0*gf(s-1,z)+gf(s-1,z+1))*wx;
  end // for z
  gf(s,np+d+1)=gf(s,np+d);
end // for s
//
//gst=sparse(gt);
//gsf=sparse(gf);
f=0.0;
endfunction

// --------------------------------------------------------------------------

function f = calc_observables(w,d)
global gst gsf np;
theta=0.0;
Zn=gf(np,1);
for z=1:d
  for s=1:np
    theta = theta + gt(s,z)*gf(np+1-s,z);
  end // for s
end // for z
f(1)=theta/w/Zn/np;
endfunction

// --------------------------------------------------------------------------

//////////////////////
//   Main program   //
//////////////////////

//////////////////////
// Global variables //
//////////////////////
global gt gf np lambda;
//////////////////
// Initial data //
//////////////////
lambda=0.16666666667; // 1/6
np=100; // number of units in the chain
d=40;
//
par0 = 0.0;
par9 = 1.5;
n=75;
pst=(par9-par0)/n;
for j=1:n+1
  par=par0+pst*(j-1);
  u=par;
  w=exp(-u);
  printf('\nj = %i   u = %f \n', j, u);
  f=calc_propagators(w,d);
  Zn=gf(np,1);
  for k=1:np
      res(k,1)=k;
      res(k,2)=gt(np,k)/Zn;
  end // for k
  // Write data into file
  s = 'pe_np=' + string(np) + '_d=' + string(d) + '_u=' + string(u) + '.dat';
  u=file('open',s,'unknown');
  write(u, res, '(1(f14.8), 32(e16.8))');
  file('close',u);
end // for j

// Plot
subplot(1,1,1); plot(res(:,1),res(:,2),'-r')
xlabel('z')
ylabel('Pe(z)')
