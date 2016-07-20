clear;
lines(0);

np=100;
q=7;
s0 = './n=' + string(np) +  '/q=' + string(q) + '/all_sig/';

j=0;

for i=175:-25:50
  j=j+1;
  sig=0.001*i
  s = s0+'SBrush_q=' + string(q) + '_' + string(i) + '.pro';
  a=fscanfMat(s);
  [m,n]=size(a);
  k=500;
  nmax=0;
  nmin=0;
  while (nmax==0)
    if(a(k-1,14) < a(k,14)) then
      nmax=1;
      kmax=k;
    end // if
    k=k-1;
  end // while
  
  while (nmin==0)
    if(a(k-1,14) > a(k,14)) then
      nmin=1;
      kmin=k;
    end // if
    k=k-1;
  end // while

  alpha=0;
  for k=kmin:m
    alpha = alpha + a(k,14)/q/sig; // nB1 (end point)
  end
//  //
  kmax
  res(j,1)=sig;
  res(j,2)=kmax;
  res(j,3)=kmin;
  res(j,4)=alpha;
  res(j,5)=(1.0-alpha)*sig;
  res(j,6)=sig*(q+1);
  res(j,7)=alpha*sig*q;
  res(j,8)=((1.0-alpha)*(q+1) + alpha)*sig;
end

// Plot
subplot(1,3,1); plot(res(:,1),res(:,2),'-or', res(:,1),res(:,3),'-og')
xlabel('sig')
ylabel('kmax, kmin')

subplot(1,3,2); plot(res(:,1),res(:,4),'-ob')
xlabel('sig')
ylabel('alpha')

subplot(1,3,3); plot(res(:,1),res(:,5),'-om')
xlabel('sig')
ylabel('(1-alpha)*sig')


/////////////////////////////////////////
// Create the name for the output file //
/////////////////////////////////////////
s0 = './';
s = s0 + 'alpha_e_SBrush_n=' + string(np) + '_q=' + string(q) + '.dat';
/////////////////////////////
// Write the data int file //
/////////////////////////////
u=file('open',s,'unknown')
write(u, res, '(f12.6,32(e16.6))')
file('close',u)

