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
  gt(s,1)=(4.0*gt(s-1,1)+gt(s-1,2))*w(1)*lambda;
  gf(s,1)=(4.0*gf(s-1,1)+gf(s-1,2))*w(1)*lambda;
  for z=2:s
    gt(s,z)=(gt(s-1,z-1)+4.0*gt(s-1,z)+gt(s-1,z+1))*w(z)*lambda;
  end // for z
  for z=2:nl-s+1
    gf(s,z)=(gf(s-1,z-1)+4.0*gf(s-1,z)+gf(s-1,z+1))*w(z)*lambda;
  end // for z
end // for s
//
f=0.0;
endfunction

// --------------------------------------------------------------------------

function f = calc_phi(w)
global gt gf nl phi;
theta=0.0;
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
//nl=150; // number of units in the chain
for nl=250:251 // number of units in the chain
q=1;
f=q+1;   // for formula M
npp=100; // for formula M
h2n=0.75; // for potention V(z)=V(h/2n)-V(z/2n), where 2n <=> npp
M=(%pi*npp/2)/acos(sqrt((f-1)/f));
//
// Calculation potention from the equation
for k=1:nl+5
    if k < h2n*2*npp then
         V(k)= -3*log(cos(h2n*%pi/2))-((3*%pi^2/8)*((h2n)^2-(h2n*2*npp/M)^2))-(-3*log(cos((%pi/2)*(k/200)))-((3*%pi^2/8)*((k/200)^2-(k/M)^2)));
    else
        V(k)=0;
    end //for if
    w(k)=exp(-V(k)); // exp(-potential) <=> stat weight
end
f=calc_propagators(w);
f=calc_phi(w);
Zn=gf(nl,1);
for k=1:nl
    res(k,1)=k;
    res(k,2)=gt(nl,k)/Zn;
    res(k,3)=phi(k);
    res(k,4)=V(k);
end // for k
//
// Write data into file
f=q+1;
s0 = '~/imc/StarBrush/Neutral/Poty/f'+ string(f) + '/h2n'+ string(h2n) +'/' ;
s = s0 + 'h2n' +string(h2n)+ '_f'+string(f)+'_z_pe_phi_m' + string(nl) + '.dat';
u=file('open',s,'unknown');
write(u, res, '(1(f14.8), 32(e16.8))');
file('close',u);
//
end // for nl
// Plot
subplot(1,3,1); plot(res(:,1),res(:,2),'-r')
xlabel('z')
ylabel('Pe(z)')

subplot(1,3,2); plot(res(:,1),res(:,3),'-b')
xlabel('z')
ylabel('phi(z)')

subplot(1,3,3); plot(res(:,1),res(:,4),'-b')
xlabel('z')
ylabel('V(z)')
