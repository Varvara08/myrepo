clear;
lines(0);

np=100;
v=0;
chi=-10;
//sig=0.1;
q=3;
//s0 = '~/varvara/Documents/IMC/n='+string(np)+'/10.2015/add/q='+string(q)+'/cs=10-1/v='+string(v)+'/';
s0 = '~/varvara/Documents/IMC/n=' + string(np) + '/chi' + string(chi)+ '/cs10-4/sig0.01/q' + string(q)+ '/v' + string(v) + '/';
j=0;

//for q=5:1:9
for sig=0.01:0.01:0.5
//  for sig=0.1
  j=j+1;
  printf('\nsig=%f', sig);
  s = s0+'SBrush_q' + string(q) +'_'+ string(j) +'.pro';
  a=fscanfMat(s);
  [m,n]=size(a);
  k=201;
 
  l=201;
  nmax=0;
  nmin=0;
  while (nmax==0)
    if(a(l-1,18) < a(l,18)) then//14/18
      nmax=1;
      lmax=l;
    end // if
    l=l-1;
  end // while
  
  while (nmin==0)
    if(a(k-1,18) > a(k,18)) then//14/18
      nmin=1;
      kmin=k;
    end // if
    k=k-1;
  end // while
  
  p=201;
  while  a(p,9)<1.0e-4 //7/9
      p=p-1;
  end

  lambda=0;
  for k=kmin:m
    lambda = lambda + a(k,18)/sig/q; // nB1 (ends point)//14/18
  end
  
    netot=0;
  for k=1:m
    netot = netot + a(k,18)/sig/q; // nB1 (ends point)//14/18
  end
//  //
  res(j,1)=sig;
  res(j,2)=lmax;
  res(j,3)=kmin; // this is h1
  res(j,4)=lambda/netot;
  res(j,5)=(1.0-lambda)*sig;
  res(j,6)=sig*(q+1);
  res(j,7)=lambda*sig*q;
  res(j,8)=((1.0-lambda)*(q+1) + lambda)*sig;
  res(j,9)=p; // this is h (total)
  h1=min(np, kmin); // this is h1 (must be <= np !!!)
  h2=p-h1; // this is h2
  res(j,10)=h1;
  res(j,11)=h2;
  res(j,12)=p/np;
  res(j,13)=h1/np;
  res(j,14)=h2/np;
  res(j,15)=lambda;
end

// Plot
subplot(1,4,1); plot(res(:,1),res(:,2),'-or', res(:,1),res(:,3),'-og')
xlabel('sig')
ylabel('kmax, kmin')

subplot(1,4,2); plot(res(:,1),res(:,4),'-og')
xlabel('sig')
ylabel('lambda')

subplot(1,4,3); plot(res(:,1),res(:,5),'-om')
xlabel('sig')
ylabel('(1-lambda)*sig')

subplot(1,4,4); plot(res(:,1),res(:,15),'-og')
xlabel('sig')
ylabel('lambda')
/////////////////////////////////////////
// Create the name for the output file //
/////////////////////////////////////////
//s0 = './';
s = s0 + 'lambda_SBrush_n' + string(np) + '_q' + string(q) + '_v' + string(v) + '_chi'+string(chi)+'_cs10-4_ends.dat';
/////////////////////////////
// Write the data int file //
/////////////////////////////
u=file('open',s,'unknown')
write(u, res, '(f12.6,32(e16.6))')
file('close',u)

