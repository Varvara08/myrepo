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

function f = calc_ends(w, np)
global gt gf g2 ends q;
Zn=0.0;
for z=1:2*np+1
    ends(z)=0.0;
    for z1=1:np
        ends(z)=ends(z)+gt(np,z1)*(gf(np+1,z1)/w(z1))^(q-1)*g2(np+1,z1,z)/w(z1);
    end
    Zn = Zn + ends(z);
end // for z
for z=1:2*np+1
    ends(z)=ends(z)/Zn;
end // for z
endfunction

// --------------------------------------------------------------------------

//////////////////////
//   Main program   //
//////////////////////

//////////////////////
// Global variables //
//////////////////////
global gt gf g2 lambda pbp ends q;
//////////////////
// Initial data //
//////////////////
lambda=0.16666666667; // 1/6
np=100; // number of units in the brush
q=2;
f=q+1;
h2n=0.75;
M=(%pi*np/2)/acos(sqrt((f-1)/f));
//
// Calculation a potention from equation  
for k=1:2*np+5
    if k < h2n*2*np then
         V(k)= -3*log(cos(h2n*%pi/2))-((3*%pi^2/8)*((h2n)^2-(h2n*2*np/M)^2))-(-3*log(cos((%pi/2)*(k/200)))-((3*%pi^2/8)*((k/200)^2-(k/M)^2)));
    else
        V(k)=0;
    end //for if
    w(k)=exp(-V(k)); // exp(-potential)
end
f=calc_gt(w,np);
f=calc_gf(w,np+1);
f=calc_g2(w,np+1);
f=calc_pbp(w,np);
f=calc_ends(w,np);
for k=1:2*np
    res(k,1)=k;
    if k < np+1 then
        res(k,2)=pbp(k);
    else
        res(k,2)=0.0;
    end
    res(k,3)=ends(k);
    res(k,4)=V(k);
    res(k,5)=w(k);

end // for k
//
// Write data into file
f=q+1;
s0 = '~/imc/StarBrush/Neutral/Poty/f'+ string(f) + '/' ;
s = s0 + 'h2n' +string(h2n)+ '_f'+string(f)+'_z_pbp_pe_np' + string(np) + '.dat';
u=file('open',s,'unknown');
write(u, res, '(1(f14.8), 32(e16.8))');
file('close',u);
//
// Plot
//subplot(1,2,1); plot(res(:,1),res(:,2),'or', res(:,1),res(:,3),'-b')
//xlabel('z')
//ylabel('pbp(z)')
//
//subplot(1,2,2); plot(res(:,1),res(:,4),'or', res(:,1),res(:,5),'-b')
//xlabel('z')
//ylabel('ends(z)')
//
// Plot
subplot(1,4,1); plot(res(:,1),res(:,2),'-r')
xlabel('z')
ylabel('Pbp(z)')

subplot(1,4,2); plot(res(:,1),res(:,3),'-b')
xlabel('z')
ylabel('end(z)')

subplot(1,4,3); plot(res(:,1),res(:,4),'-b')
xlabel('z')
ylabel('V(z)')

//subplot(1,4,4); plot(res(:,1),res(:,5),'-b')
//xlabel('z')
//ylabel('w(z)')
//

