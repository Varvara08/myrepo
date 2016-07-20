clear;
clear gt gf;
lines(0);

// --------------------------------------------------------------------------

function f = calc_gt(w, np)
global gt lambda;
//
// Array initialization with zeros
for s=1:np
  for z=1:np+1
    gt(s,z)=0.0;
  end // for z
end // for s
//
[m,n]=size(w);
// Initial condition
gt(1,1)=w(1);

// Start recurrence
for s=2:np
  // Boundary condition (2 lines) 
  gt(s,1)=(4.0*gt(s-1,1)+gt(s-1,2))*w(1)*lambda;
  for z=2:s
    gt(s,z)=(gt(s-1,z-1)+4.0*gt(s-1,z)+gt(s-1,z+1))*w(z)*lambda;
  end // for z
end // for s
//
f=0.0;
endfunction

// --------------------------------------------------------------------------

function f = calc_gf(w, np)
global gf lambda;
//
// Array initialization with zeros
for s=1:np
  for z=1:2*np+1
    gf(s,z)=0.0;
  end // for z
end // for s
//
[m,n]=size(w);
// Initial condition
for z=1:2*np+1
  gf(1,z)=w(z);
end // for z

// Start recurrence
for s=2:np
  // Boundary condition (2 lines) 
  gf(s,1)=(4.0*gf(s-1,1)+gf(s-1,2))*w(1);
  for z=2:2*np-s+1
    gf(s,z)=(gf(s-1,z-1)+4.0*gf(s-1,z)+gf(s-1,z+1))*w(z);
  end // for z
end // for s
//
f=0.0;
endfunction

// --------------------------------------------------------------------------

function f = calc_g2(w, np)
global g2 lambda;
//
// Array initialization with zeros
for s=1:np
  for z1=1:np
    for z=1:2*np+1
        g2(s,z1,z)=0.0;
    end
  end // for z
end // for s
//
// Initial condition
for z1=1:np
    g2(1,z1,z1)=w(z1);
end
// Start recurrence
for s=2:np
    for z1=1:np
        // Boundary condition
        g2(s,z1,1)=(4.0*g2(s-1,z1,1)+g2(s-1,z1,2))*w(1)*lambda;
        for z=2:s+z1-1
            g2(s,z1,z)=(g2(s-1,z1,z-1)+4.0*g2(s-1,z1,z)+g2(s-1,z1,z+1))*w(z)*lambda;
        end // for z
    end // for z1
end // for s
//
f=0.0;
endfunction

// --------------------------------------------------------------------------

function f = calc_pbp(w, np)
global gt gf pbp q;
Zn=0.0;
for z=1:np
    pbp(z)=gt(np,z)*(gf(np+1,z)/w(z))^q;
    Zn = Zn + pbp(z);
end // for z
for z=1:np
    pbp(z)=pbp(z)/Zn;
end // for z
endfunction

// --------------------------------------------------------------------------

function f = calc_pbp2(w, np)
global gt gf g2 pbp2 q;
Zn=0.0;
for z=1:2*np+1
    pbp2(z)=0.0;
    for z1=1:np
        pbp2(z)=pbp2(z)+gt(np,z1)*(gf(np+1,z1)/w(z1))^(q-1)*g2(np+1,z1,z)/w(z1);
    end
    Zn = Zn + pbp2(z);
end // for z
for z=1:2*np+1
    pbp2(z)=pbp2(z)/Zn;
end // for z
endfunction

// --------------------------------------------------------------------------

//////////////////////
//   Main program   //
//////////////////////

//////////////////////
// Global variables //
//////////////////////
global gt gf g2 lambda pbp pbp2 q;
//////////////////
// Initial data //
//////////////////
lambda=0.16666666667; // 1/6
np=100; // number of units in the brush
q=2;
f=q+1;
//
// Открываем --> считаем потенциал щётки для различного слоя одной сигмы. к-слой
s0 = '~/imc/StarBrush/Neutral/Poty/f'+ string(f) + '/' ;
s = s0 +'CSBrush_m240_300.pro';
a=fscanfMat(s);
[m,n]=size(a);
tbp=0;
tbp2=0;
for k=1:m
    w(k)=1-a(k,7); // exp(-potential)
    tbp = tbp + a(k,13);
    tbp2 = tbp2 + a(k,15);
end
f=calc_gt(w,np);
f=calc_gf(w,np+1);
f=calc_g2(w,np+1);
f=calc_pbp(w,np);
f=calc_pbp2(w,np);
for k=1:2*np
    res(k,1)=k;
    if k < np+1 then
        res(k,2)=pbp(k);
    else
        res(k,2)=0.0;
    end
    res(k,3)=a(k,13)/tbp;
    res(k,4)=pbp2(k);
    res(k,5)=a(k,15)/tbp2;
end // for k
//
// Write data into file
f=q+1;
s = s0 + 'Pf'+string(f)+'_z_pbp_pe_np' + string(np) + '.dat';
u=file('open',s,'unknown');
write(u, res, '(1(f14.8), 32(e16.8))');
file('close',u);
//
// Plot
subplot(1,2,1); plot(res(:,1),res(:,2),'or', res(:,1),res(:,3),'-b')
xlabel('z')
ylabel('pbp(z)')

subplot(1,2,2); plot(res(:,1),res(:,4),'or', res(:,1),res(:,5),'-b')
xlabel('z')
ylabel('pbp2(z)')

