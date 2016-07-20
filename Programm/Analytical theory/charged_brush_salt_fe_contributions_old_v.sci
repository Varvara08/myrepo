clear;
lines(0);

// --------------------------------------------------------------------------

function f = fe(y) 
global q sig chi alpha cs
//
x = y(1);
lambda = y(2);


g = -log((2+cosh(x))/3) - log((2+cosh(x/q))/3)*q;

dz1=sinh(x)/(2+cosh(x)); // степень растяжения стебля
dz2=sinh(x/q)/(2+cosh(x/q)); // степень растяжения свободных ветвей 

phi1 = ((q+1)-lambda*q)*sig/dz1; // объемная доля полимера в нижнем слое
phi2 = lambda*q*sig/dz2; // объемная доля полимера в верхнем слое

v1 = (1-phi1)/phi1*log(1-phi1) + chi*(1-phi1) + (1-(1+(alpha*phi1/cs)^2)^0.5)/phi1*cs + alpha*asinh(alpha*phi1/cs);
v2 = (1-phi2)/phi2*log(1-phi2) + chi*(1-phi2) + (1-(1+(alpha*phi2/cs)^2)^0.5)/phi2*cs + alpha*asinh(alpha*phi2/cs);

f = lambda*(g + dz1*x + dz2*x); // конформационная свободная энергия растянутых звезд
f = f + (lambda+(q+1)*(1-lambda))*v1+lambda*q*v2; // + вклады  в каждом слое от трансляционной энтропии растворенных молекул и взаимодействия полимер-растворитель и полимер-полимер, осмотическое давление из-за разницы концентраций ионов внутри и вне щетки,а также энтропия, связанная с распределением заряженных и незаряженных сегментов в пэ цепи 

endfunction
// --------------------------------------------------------------------------
function f = fobs(y) //h and sig and F 
global q sig chi alpha cs
//
x = y(1);
lambda = y(2);

g = -log((2+cosh(x))/3) - log((2+cosh(x/q))/3)*q;

dz1=sinh(x)/(2+cosh(x)); // степень растяжения стебля
dz2=sinh(x/q)/(2+cosh(x/q)); // степень растяжения свободных ветвей 

phi1 = ((q+1)-lambda*q)*sig/dz1; // объемная доля полимера в нижнем слое
phi2 = lambda*q*sig/dz2; // объемная доля полимера в верхнем слое

v1 = (1-phi1)/phi1*log(1-phi1) + chi*(1-phi1) + (1-(1+(alpha*phi1/cs)^2)^0.5)/phi1*cs + alpha*asinh(alpha*phi1/cs);
v2 = (1-phi2)/phi2*log(1-phi2) + chi*(1-phi2) + (1-(1+(alpha*phi2/cs)^2)^0.5)/phi2*cs + alpha*asinh(alpha*phi2/cs);

ff = fe(y);

f(1)=dz1;
f(2)=dz2;
f(3)=(dz1+dz2)/2;
f(4)=phi1;
f(5)=phi2;
f(6)=(q+1)*sig/(dz1+dz2);//средняя phi
f(7)=ff;
f(8)=lambda*(g + dz1*x + dz2*x); // конформационная свободная энергия растянутых звезд
f(9)=(lambda+(q+1)*(1-lambda))*(1-phi1)/phi1*log(1-phi1);//вклад от трансляционной энтропии растворителя в первом слое
f(10)=lambda*q*(1-phi2)/phi2*log(1-phi2);//вклад от трансляционной энтропии растворителя во втором слое
f(11)=(lambda+(q+1)*(1-lambda))*chi*(1-phi1);//взаимодействие полимер-полимер и полимер-растворитель в первом слое
f(12)=lambda*q*chi*(1-phi2);//взаимодействие полимер-полимер и полимер-растворитель во втором слое
f(13)=(lambda+(q+1)*(1-lambda))*(1-(1+(alpha*phi1/cs)^2)^0.5)/phi1*cs;//вклад от энтропии,связанной с распределением заряженных и незаряженных сегментов пэ цепи в первом слое
f(14)=lambda*q*(1-(1+(alpha*phi2/cs)^2)^0.5)/phi2*cs;//вклад от энтропии,связанной с распределением заряженных и незаряженных сегментов пэ цепи во втором слое
f(15)=(lambda+(q+1)*(1-lambda))*alpha*asinh(alpha*phi1/cs);//осмотическое давление, обусловленное различием в ионных концентрациях между внутренними и внешними частями щетки в первом слое
f(16)=lambda*q*alpha*asinh(alpha*phi2/cs);//осмотическое давление, обусловленное различием в ионных концентрациях между внутренними и внешними частями щетки во втором слое

endfunction
// --------------------------------------------------------------------------

function f = dfe(y) 
global q sig chi alpha cs
//
x = y(1);
lambda = y(2);


g = -log((2+cosh(x))/3) - log((2+cosh(x/q))/3)*q;

dz1=sinh(x)/(2+cosh(x)); // степень растяжения стебля
dz2=sinh(x/q)/(2+cosh(x/q)); // степень растяжения свободных ветвей 

phi1 = ((q+1)-lambda*q)*sig/dz1; // объемная доля полимера в нижнем слое
phi2 = lambda*q*sig/dz2; // объемная доля полимера в верхнем слое

v1 = (1-phi1)/phi1*log(1-phi1) + chi*(1-phi1) + (1-(1+(alpha*phi1/cs)^2)^0.5)/phi1*cs + alpha*asinh(alpha*phi1/cs);
v2 = (1-phi2)/phi2*log(1-phi2) + chi*(1-phi2) + (1-(1+(alpha*phi2/cs)^2)^0.5)/phi2*cs + alpha*asinh(alpha*phi2/cs);

dfedlam = g + dz1*x + dz2*x;
dfedlam = dfedlam - q*v1 + q*v2;

dfedphi1 =( -log(1-phi1)/phi1/phi1 - 1/phi1 - chi + alpha^2/cs/(1+(alpha*phi1/cs)^2)^0.5 + cs/(1+(alpha*phi1/cs)^2)^0.5/phi1^2 - cs/phi1^2 )*(lambda+(q+1)*(1-lambda));

dfedphi2 = (-log(1-phi2)/phi2/phi2 - 1/phi2 - chi + alpha^2/cs/(1+(alpha*phi2/cs)^2)^0.5 + cs/(1+(alpha*phi2/cs)^2)^0.5/phi2^2 - cs/phi2^2)*lambda*q;


dphi1dlam = -q*sig/dz1;
dphi2dlam = q*sig/dz2;

f(1) = dfedlam+dfedphi1*dphi1dlam+dfedphi2*dphi2dlam;

dfeddz1 = x*lambda;
dfeddz2 = x*lambda;

ddz1dx = (2*cosh(x)+1.0)/(2+cosh(x))^2;
ddz2dx = (2*cosh(x/q)+1.0)/q/(2+cosh(x/q))^2;
dphi1ddz1 = -phi1/dz1;
dphi2ddz2 = -phi2/dz2;

dfedx = lambda*(- sinh(x)/(2+cosh(x)) - sinh(x/q)/(2+cosh(x/q)) + dz1 + dz2);

f(2) = dfedx + dfeddz1*ddz1dx + dfeddz2*ddz2dx
f(2) = f(2) + dfedphi1*dphi1ddz1*ddz1dx + dfedphi2*dphi2ddz2*ddz2dx;



endfunction

//////////////////
// Main program //
//////////////////

//////////////////////
// Global variables //
//////////////////////
global q sig chi alpha cs

q = 2 ;
chi=0;
alpha=0;
cs=10^(-4);
y0 = [4; 0.3]
par0 = 0.2;
par9 = 0.4;
n=100;

//y0 = [6.21; 0.68]//
//
//par0 = 0.6;
//par9 = 0.26;
//n=200;

pst=(par9-par0)/n;
for k=1:n+1
  par=par0+pst*(k-1);
  sig=par;
    printf("\nsig=%f", sig);
  res(k,1)=par;
  [ymin,fval,info]=fsolve(y0,dfe);
  y0=ymin;                              //golden rule
  res(k,2)=ymin(1);
  res(k,3)=ymin(2);
  res(k,4)=fval(1);
  res(k,5)=fval(2);
  [J,H]=derivative(fe,ymin,H_form='hypermat');
  res(k,6)=J(1);
  res(k,7)=J(2);
  res(k,8)=det(H);
  obs=fobs(ymin)
  res(k,9)=obs(1);
  res(k,10)=obs(2);
  res(k,11)=obs(3);
  res(k,12)=obs(4);
  res(k,13)=obs(5);
  res(k,14)=obs(6);
  res(k,15)=obs(7);
  res(k,16)=obs(8);
  res(k,17)=obs(9);
  res(k,18)=obs(10);
  res(k,19)=obs(11);
  res(k,20)=obs(12);
  res(k,21)=obs(13);
  res(k,22)=obs(14);
  res(k,23)=obs(15);
  res(k,24)=obs(16);
  res(k,25)=ymin(2)/(1-ymin(2));
end


s0 = '~/varvara/Documents/IMC/Analytical theory/Charged/before_ex/';
sout = 'lambda_vs_sig_q=' + string(q) + '_v'+ string(alpha)+'_v1.dat';
u=file('open',sout,'unknown');
write(u, res, '(f12.6,32(e16.6))');
file('close',u);

// Plot
subplot(4,4,1); plot(res(:,1),res(:,2),'-r')
xlabel('sig')
ylabel('force_1')
//
subplot(4,4,2); plot(res(:,1),res(:,3),'-r')
xlabel('sig')
ylabel('lambda')
//
subplot(4,4,3); plot(res(:,1),res(:,4),'-r',res(:,1),res(:,5),'-g')
xlabel('sig')
ylabel('fval')
//
// Plot
subplot(4,4,4); plot(res(:,1),res(:,16),'-r')
xlabel('sig')
ylabel('force_conf.star')
//
// Plot
subplot(4,4,5); plot(res(:,1),res(:,17),'-r')
xlabel('sig')
ylabel('force_sol_1')
//
// Plot
subplot(4,4,6); plot(res(:,1),res(:,18),'-r')
xlabel('sig')
ylabel('force_sol_2')
//
// Plot
subplot(4,4,7); plot(res(:,1),res(:,19),'-r')
xlabel('sig')
ylabel('force_int_1')
//
// Plot
subplot(4,4,8); plot(res(:,1),res(:,20),'-r')
xlabel('sig')
ylabel('force_int_2')
//
// Plot
subplot(4,4,9); plot(res(:,1),res(:,21),'-r')
xlabel('sig')
ylabel('force_ch_1')
//
// Plot
subplot(4,4,10); plot(res(:,1),res(:,22),'-r')
xlabel('sig')
ylabel('force_ch_2')
//
// Plot
subplot(4,4,11); plot(res(:,1),res(:,23),'-r')
xlabel('sig')
ylabel('force_os.pr_1')
//
// Plot
subplot(4,4,12); plot(res(:,1),res(:,24),'-r')
xlabel('sig')
ylabel('force_os.pr_2')
