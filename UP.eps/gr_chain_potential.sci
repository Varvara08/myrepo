clear;
clear gt gf;
lines(0);

function f = calc_propagators(w)
global gt gf nl lambda;
//
// Array initialization with zeros
// gt
for s=1:nl
  for z=1:nl+1
    gt(s,z)=0.0;
    gf(s,z)=0.0;
  end // for z
end // for s
// gf
//
[m,n]=size(w);
// Initial condition
gt(1,1)=w(1);
for z=1:nl
  gf(1,z)=w(z);
end // for z

// Start recurrence
for s=2:nl
  gt(s,1)=(4.0*gt(s-1,1)+gt(s-1,2))*w(1);
  gf(s,1)=(4.0*gf(s-1,1)+gf(s-1,2))*w(1);
  for z=2:s
    gt(s,z)=(gt(s-1,z-1)+4.0*gt(s-1,z)+gt(s-1,z+1))*w(z);
  end // for z
  for z=2:nl-s+1
    gf(s,z)=(gf(s-1,z-1)+4.0*gf(s-1,z)+gf(s-1,z+1))*w(z);
  end // for z
end // for s
//
f=0.0;
endfunction

// --------------------------------------------------------------------------

function f = calc_phi(w)
global gt gf nl phi;
Zn=gf(nl,1);
for z=1:nl
    phi(z)=0.0;
    for s=1:nl
        phi(z) = phi(z) + gt(s,z)*gf(nl+1-s,z)/w(z);
    end // for s
end // for z
for z=1:nl
    phi(z)=phi(z)/Zn;
end // for z
endfunction

// --------------------------------------------------------------------------

//////////////////////
//   Main program   //
//////////////////////

//////////////////////
// Global variables //
//////////////////////
global gt gf nl lambda phi;
//////////////////
// Initial data //
//////////////////
lambda=0.16666666667; // 1/6
nl=260; // number of units in the chain
np=240; // for opening file
q=2;
f=q+1;
//sig=0.1*100;
//
// Открываем --> считаем потенциал щётки для различного слоя одной сигмы. к-слой
s0 = '~/imc/StarBrush/Neutral/Poty/f'+ string(f) + '/' ;
s = s0 +'CSBrush_m' +string(np)+'_100.pro';
a=fscanfMat(s);
[m,n]=size(a);
for k=1:m
    w(k)=1-a(k,7); // exp(-potential) <=> statistical weight
end
f=calc_propagators(w);
f=calc_phi(w);
Zn=gf(nl,1);
for k=1:nl
    res(k,1)=k;
    res(k,2)=phi(k);
    res(k,3)=gt(nl,k)/Zn;
end // for k
//
// Write data into file
f=q+1;
s = s0 + 'Pf'+string(f)+'_z_phi_pe_m' + string(nl) + '.dat';
u=file('open',s,'unknown');
write(u, res, '(1(f14.8), 32(e16.8))');
file('close',u);
//
// Plot
subplot(1,2,1); plot(res(:,1),res(:,2),'-r')
xlabel('z')
ylabel('Pe(z)')

subplot(1,2,2); plot(res(:,1),res(:,3),'-b')
xlabel('z')
ylabel('phi(z)')
