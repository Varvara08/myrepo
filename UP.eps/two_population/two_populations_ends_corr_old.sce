clear;
lines(0);

np=100;
v=0;
chi=-10;
//sig=0.1;
q=3;
//s0 = '~/varvara/Documents/IMC/n='+string(np)+'/10.2015/add/q='+string(q)+ '/cs=10-1/v='+string(v)+'/';
//s0 = '~/varvara/Documents/IMC/n='+string(np)+'/10.2015/add/sig=0.1_cs=10-1/q='+string(q)+ '/v='+string(v)+'/';
//s0 = '~/varvara/Documents/IMC/n='+string(np)+'/0/cs=10-1/q='+string(q)+ '/v='+string(v)+'/';
//s0 = '~/varvara/Documents/IMC/n='+string(np)+'/cs=10-4/old/q='+string(q)+ '/phi_s=10-4/v='+string(v)+'/';
s0 = '~/varvara/Documents/IMC/n=' + string(np) + '/chi' + string(chi)+ '/cs10-4/sig0.01/q' + string(q)+ '/v' + string(v) + '/';
j=0;

//for q=5:1:9
for sig=0.01:0.01:0.5
  //for sig=0.1
  j=j+1;
  printf('\nsig=%f', sig);
 s = s0+'SBrush_q' + string(q) + '_'+string(j)+ '.pro';
  //s = s0+'SBrush_q=' + string(q) + '.pro';
  a=fscanfMat(s);
  [m,n]=size(a);
  k=201;
  p=201;

  nmax=0;
  nmin=0;
  while (nmax==0)&(k>1)
    if(a(k-1,18) < a(k,18)) then
      nmax=1;
      kmax=k;
    end // if
    k=k-1;
  end // while
  
  k=k+1;
  
  while (nmin==0)&(k>1)
    if(a(k-1,18) > a(k,18)) then
      nmin=1;
      kmin=k;
    end // if
    k=k-1;
  end // while
  
  if k == 1 then
      kmin=1;
  end
  
  while  a(p,9)<1.0e-4 
      p=p-1;
  end

  lambda=0;
  for k=kmin:m
    lambda = lambda + a(k,18)/sig/q; // nB1 (ends point)
  end
  
    netot=0;
  for k=1:m
    netot = netot + a(k,18)/sig/q; // nB1 (ends point)
  end
//  //
  res(j,1)=sig;
  res(j,2)=kmax;
  res(j,3)=kmin; // this is h1
  res(j,4)=lambda/netot;
  res(j,5)=(1.0-lambda)*sig;
  res(j,6)=sig*(q+1);
  res(j,7)=lambda*sig*q;
  res(j,8)=((1.0-lambda)*(q+1) + lambda)*sig;
  res(j,9)=p; // this is h (total)
  if kmin>1 then
      //h1=min(np, kmin); // this is h1 (must be <= np !!!)
      
       h1=kmin; // this is h1 (must be <= np !!!)
       h2=p-h1; // this is h2
  else
      h1=p; // this is h1 (must be <= np !!!)
      h2=0;
  end
  
  res(j,10)=h1;
  res(j,11)=h2;
  res(j,12)=p/np;
  res(j,13)=h1/np;
  res(j,14)=h2/np;
end

// Plot
subplot(2,2,1); plot(res(:,1),res(:,2),'-or', res(:,1),res(:,3),'-og')
xlabel('sig')
ylabel('kmax, kmin')

subplot(2,2,2); plot(res(:,1),res(:,4),'-og')
xlabel('sig')
ylabel('lambda')

subplot(2,2,3); plot(res(:,1),res(:,9),'-or',res(:,1),res(:,10),'-og',res(:,1),res(:,11),'-om')
xlabel('sig')
ylabel('H,H_1,H_2')

subplot(2,2,4); plot(res(:,9),res(:,10),'-og',res(:,9),res(:,11),'-om')
xlabel('H')
ylabel('H_1,H_2')



/////////////////////////////////////////
// Create the name for the output file //
/////////////////////////////////////////
//s0 = './';
s = s0 + 'lambda_SBrush_n' + string(np) +'_q' +string(q)+'_v' + string(v) + '_chi'+string(chi)+'_cs10-4_ends.dat';
/////////////////////////////
// Write the data int file //
/////////////////////////////
u=file('open',s,'unknown')
write(u, res, '(f12.6,32(e16.6))')
file('close',u)

