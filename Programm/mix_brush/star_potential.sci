clear;
clear gt gf;
lines(0);

function f = calc_propagators(w, np)
global gt gf lambda;
//
// Array initialization with zeros
// gt
for s=1:np
  for z=1:np+1
    gt(s,z)=0.0;
    gf(s,z)=0.0;
  end // for z
end // for s
// gf
//
[m,n]=size(w);
// Initial condition
gt(1,1)=w(1);
for z=1:np
  gf(1,z)=w(z);
end // for z

// Start recurrence
for s=2:np
  gt(s,1)=(4.0*gt(s-1,1)+gt(s-1,2))*w(1);
  gf(s,1)=(4.0*gf(s-1,1)+gf(s-1,2))*w(1);
  for z=2:s
    gt(s,z)=(gt(s-1,z-1)+4.0*gt(s-1,z)+gt(s-1,z+1))*w(z);
  end // for z
  for z=2:np-s+1
    gf(s,z)=(gf(s-1,z-1)+4.0*gf(s-1,z)+gf(s-1,z+1))*w(z);
  end // for z
end // for s
//
f=0.0;
endfunction

// --------------------------------------------------------------------------

function f = calc_pbp(w, np)
global gt gf pbp q;
Zn=0.0;
for z=1:np
    pbp(z)=gt(np,z)*(gf(np,z)/w(z))^q;
    Zn = Zn + pbp(z);
end // for z
for z=1:np
    pbp(z)=pbp(z)/Zn;
end // for z
endfunction

// --------------------------------------------------------------------------

//////////////////////
//   Main program   //
//////////////////////

//////////////////////
// Global variables //
//////////////////////
global gt gf lambda pbp q;
//////////////////
// Initial data //
//////////////////
lambda=0.16666666667; // 1/6
np=100; // number of units in the chain
q=1;
//
// Открываем --> считаем потенциал щётки для различного слоя одной сигмы. к-слой
s = 'SBrush_q=1_300.pro';
a=fscanfMat(s);
[m,n]=size(a);
tbp=0;
for k=1:m
    w(k)=1-a(k,7); // exp(-potential)
    tbp = tbp + a(k,12);
end
f=calc_propagators(w,m);
f=calc_pbp(w,np);
for k=1:np
    res(k,1)=k;
    res(k,2)=pbp(k);
    res(k,3)=a(k,12)/tbp;
end // for k
//
// Write data into file
s = 'pbp_np=' + string(np) + '.dat';
u=file('open',s,'unknown');
write(u, res, '(1(f14.8), 32(e16.8))');
file('close',u);
//
// Plot
subplot(1,1,1); plot(res(:,1),res(:,2),'-r', res(:,1),res(:,3),'-b')
xlabel('z')
ylabel('pbp(z)')
