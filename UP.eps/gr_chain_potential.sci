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
Zn=gf(nl,1);
for z=1:nl
    phi(z)=0.0;
    for s=1:nl
        phi(z) = phi(z) + gt(s,z)*gf(nl+1-s,z)/w(z);
    end // for s
end // for z
//for z=1:nl
//    phi(z)=phi(z)/Zn;
//end // for z
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
////nl=260; // number of units in the chain
///////My cycle for parameter///
for nl=250
////////////////////////////////
//np=240;  // for opening file
q=1;
f=q+1;   // for output.file name
h2n=0.25 // for output.file name
M = 435  // for output.file name this M in дроби
//sig=0.1*100;
//
// Открываем --> считаем потенциал щётки для различного слоя одной сигмы. к-слой
//s0 = '~/imc/StarBrush/Neutral/Poty/f'+ string(f) + '/' ; // (1)
s0 = '~/imc/StarBrush/Neutral/Poty/simpleFiles/' ; // (brush)
//s = s0 +'CSBrush_m' +string(np)+'_100.pro';
//s = s0 +'output_h2n'+string(h2n)+'_f'+string(f)+'_sig950_z_pbp_pe_V1_np100_'+string(M)+'.dat'; (1)
//s = s0 +'sfbox_f3_950_man_munis_cutted_h2n0.25_'+string(M)+'.dat'; // (2)
//s = s0 +'sfbox_f3_95.dat'; // (2)
s = s0 +'my_phi.dat'; // (2)
a=fscanfMat(s);
[m,n]=size(a);
for k=1:nl
    if k <= m then
        V(k)=-log(1-a(k,1)) ////// a(k,4) --- (1) ; a(k,2) --- (2)
    else 
        V(k)=0;
    end // for if
//w(k)=1-a(k,7); // exp(-potential) <=> statistical weight
w(k)=exp(-V(k)); // exp(-potential) <=> statistical weight in output.file

end
f=calc_propagators(w);
f=calc_phi(w);
Zn=gf(nl,1);
for k=1:nl
    res(k,1)=k;
    res(k,2)=phi(k);
    res(k,3)=gt(nl,k)/Zn;
    res(k,4)=V(k);
end // for k
//
// Write data into file
f=q+1; // for output name
s0 = '~/imc/StarBrush/Neutral/Poty/f'+ string(f) + '/' ; // (2)
//s = s0 +'/h2n'+string(h2n)+'/'+'h2n'+string(h2n)+'_f'+string(f)+'_z_phi_pe_m' + string(nl) + '.dat'; (1)
//s = s0 +'/h2n'+string(h2n)+'/M'+string(M)+'/'+'h2n'+string(h2n)+'_f'+string(f)+'8_z_phi_pe_m' + string(nl) + '.dat'; // (2)
s = s0 +'mytest_z_phi_pe_V.dat'; // (2)
u=file('open',s,'unknown');
write(u, res, '(1(f14.8), 32(e16.8))');
file('close',u);
//
// Plot
subplot(1,3,1); plot(res(:,1),res(:,2),'-r')
xlabel('z')
ylabel('phi(z)')

subplot(1,3,2); plot(res(:,1),res(:,3),'-b')
xlabel('z')
ylabel('Pe(z)')

subplot(1,3,3); plot(res(:,1),res(:,4),'-g')
xlabel('z')
ylabel('V(z)')

/////////////////End my cycle///
end // for for
////////////////////////////////
