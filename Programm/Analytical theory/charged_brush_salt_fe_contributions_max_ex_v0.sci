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
dz3 = dz1/2.0; // растяжение длинного пути
dz4 = dz2/2.0;

s = sqrt(1.0 + 3.0*dz3^2) - 1.0;


phi1 = ((q+1)-lambda*q)*sig/dz1; // объемная доля полимера в нижнем слое
phi2 = lambda*q*sig/dz2; // объемная доля полимера в верхнем слое

v1 = (1-phi1-alpha*phi1)/phi1*log(1-phi1-alpha*phi1) + chi*(1-phi1) + (1-(1+(alpha*phi1/cs)^2)^0.5)/phi1*cs + alpha*asinh(alpha*phi1/cs);
v2 = (1-phi2-alpha*phi2)/phi2*log(1-phi2-alpha*phi2) + chi*(1-phi2) + (1-(1+(alpha*phi2/cs)^2)^0.5)/phi2*cs + alpha*asinh(alpha*phi2/cs);

f = lambda*(g + dz1*x + dz2*x)/(q+1); // конформационная свободная энергия растянутых звезд
// ?!??! делить на q+1

f = f +(1.0-lambda)/(q+1) *((1+2*s2+3*dz4)/6 *log(1+2*s2+3*dz4) + (1+2*s2-3*dz4)/6 *log(1+2*s2-3*dz4) + 2/3 *(1-s2)*log(1-s2));///
f = f +(1.0-lambda)/(q+1) *((1+2*s1+3*dz3)/6 *log(1+2*s1+3*dz3) + (1+2*s1-3*dz3)/6 *log(1+2*s1-3*dz3) + 2/3 *(1-s1)*log(1-s1));// + вклад от растяжения длинного пути звезд нижней популяции

f = f + (lambda+(q+1)*(1-lambda))/(q+1)*v1+lambda*q/(q+1)*v2; // + вклады  в каждом слое от трансляционной энтропии растворенных молекул и взаимодействия полимер-растворитель и полимер-полимер, осмотическое давление из-за разницы концентраций ионов внутри и вне щетки,а также энтропия, связанная с распределением заряженных и незаряженных сегментов в пэ цепи 
//?!??! делить на q+1

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
dz3 = dz1/2.0; //  растяжение длинного пути
dz4 = dz2/2.0;


s1 = sqrt(1.0 + 3.0*dz3^2) - 1.0;
s2 = sqrt(1.0 + 3.0*dz4^2) - 1.0;////

phi1 = ((q+1)-lambda*q)*sig/dz1; // объемная доля полимера в нижнем слое
phi2 = lambda*q*sig/dz2; // объемная доля полимера в верхнем слое

v1 = (1-phi1-alpha*phi1)/phi1*log(1-phi1-alpha*phi1) + chi*(1-phi1) + (1-(1+(alpha*phi1/cs)^2)^0.5)/phi1*cs + alpha*asinh(alpha*phi1/cs);
v2 = (1-phi2-alpha*phi2)/phi2*log(1-phi2-alpha*phi2) + chi*(1-phi2) + (1-(1+(alpha*phi2/cs)^2)^0.5)/phi2*cs + alpha*asinh(alpha*phi2/cs);

ff = fe(y);

f(1)=dz1;
f(2)=dz2;
f(3)=(dz1+dz2)/2;
f(4)=phi1;
f(5)=phi2;
f(6)=(q+1)*sig/(dz1+dz2);//средняя phi
f(7)=ff;

// ?!??! делить на q+1

f(8)=lambda*(g + dz1*x + dz2*x)/(q+1); // конформационная свободная энергия растянутых звезд
f(9)=(lambda+(q+1)*(1-lambda))/(q+1)*(1-phi1-alpha*phi1)/phi1*log(1-phi1-alpha*phi1);//вклад от трансляционной энтропии растворителя в первом слое
f(10)=lambda*q/(q+1)*(1-phi2-alpha*phi2)/phi2*log(1-phi2-alpha*phi2);//вклад от трансляционной энтропии растворителя во втором слое
f(11)=(lambda+(q+1)*(1-lambda))/(q+1)*chi*(1-phi1);//взаимодействие полимер-полимер и полимер-растворитель в первом слое
f(12)=lambda*q/(q+1)*chi*(1-phi2);//взаимодействие полимер-полимер и полимер-растворитель во втором слое
f(13)=(lambda+(q+1)*(1-lambda))/(q+1)*(1-(1+(alpha*phi1/cs)^2)^0.5)/phi1*cs;//вклад от энтропии,связанной с распределением заряженных и незаряженных сегментов пэ цепи в первом слое
f(14)=lambda*q/(q+1)*(1-(1+(alpha*phi2/cs)^2)^0.5)/phi2*cs;//вклад от энтропии,связанной с распределением заряженных и незаряженных сегментов пэ цепи во втором слое
f(15)=(lambda+(q+1)*(1-lambda))/(q+1)*alpha*asinh(alpha*phi1/cs);//осмотическое давление, обусловленное различием в ионных концентрациях между внутренними и внешними частями щетки в первом слое
f(16)=lambda*q/(q+1)*alpha*asinh(alpha*phi2/cs);//осмотическое давление, обусловленное различием в ионных концентрациях между внутренними и внешними частями щетки во втором слое

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
dz3 = dz1/2.0; //  растяжение длинного пути
dz4 = dz2/2.0;


s1 = sqrt(1.0 + 3.0*dz3^2) - 1.0;
s2 = sqrt(1.0 + 3.0*dz4^2) - 1.0;////

phi1 = ((q+1)-lambda*q)*sig/dz1; // объемная доля полимера в нижнем слое
phi2 = lambda*q*sig/dz2; // объемная доля полимера в верхнем слое

v1 = (1-phi1-alpha*phi1)/phi1*log(1-phi1-alpha*phi1) + chi*(1-phi1) + (1-(1+(alpha*phi1/cs)^2)^0.5)/phi1*cs + alpha*asinh(alpha*phi1/cs);
v2 = (1-phi2-alpha*phi2)/phi2*log(1-phi2-alpha*phi2) + chi*(1-phi2) + (1-(1+(alpha*phi2/cs)^2)^0.5)/phi2*cs + alpha*asinh(alpha*phi2/cs);

//dz

ddz1dx = (2*cosh(x)+1.0)/(2+cosh(x))^2;
ddz2dx = (2*cosh(x/q)+1.0)/q/(2+cosh(x/q))^2;
ddz3dx = ddz1dx/2.0;
ddz4dx = ddz2dx/2.0;////

//phi
dphi1dlam = -q*sig/dz1;
dphi2dlam = q*sig/dz2;

dphi1ddz1 = -phi1/dz1;
dphi2ddz2 = -phi2/dz2;

//s
ds1ddz3 = 3*dz3/sqrt(1.0+3.0*dz3^2);
ds2ddz4 = 3*dz4/sqrt(1.0+3.0*dz4^2);////

//f(1)
dfedx = lambda*x/(q+1) * (ddz1dx + ddz2dx);
dfedx = dfedx + (1-lambda)/(q+1) * ((1/3 *ds1ddz3+0.5)*log(1+2*s1+3*dz3) + (ds1ddz3/3 -0.5)*log(1+2*s1-3*dz3) - 2/3 *ds1ddz3*log(1-s1))*ddz3dx;
dfedx = dfedx + (1-lambda)/(q+1) * ((1/3 *ds2ddz4+0.5)*log(1+2*s2+3*dz4) + (ds2ddz4/3 -0.5)*log(1+2*s2-3*dz4) - 2/3 *ds2ddz4*log(1-s2))*ddz4dx;///
dfedx = dfedx + (q*(1-lambda)+1)/(q+1) * dphi1ddz1*ddz1dx*(-log(1-phi1-alpha*phi1)/phi1^2 - (1+alpha)/phi1 - chi - cs/phi1^2 + cs/phi1^2 *(1+(alpha*phi1/cs)^2)^0.5);
dfedx = dfedx + q*lambda/(q+1) * dphi2ddz2*ddz2dx*(-log(1-phi2-alpha*phi2)/phi2^2 - (1+alpha)/phi2 - chi - cs/phi2^2 + cs/phi2^2 *(1+(alpha*phi2/cs)^2)^0.5);

//f(2)
fa = (1+2*s1+3*dz3)/6 *log(1+2*s1+3*dz3) + (1+2*s1-3*dz3)/6 *log(1+2*s1-3*dz3) + 2/3 *(1-s1)*log(1-s1);
fa = fa + (1+2*s2+3*dz4)/6 *log(1+2*s2+3*dz4) + (1+2*s2-3*dz4)/6 *log(1+2*s2-3*dz4) + 2/3 *(1-s2)*log(1-s2);////

dfedlambda = (g + dz1*x + dz2*x)/(q+1) - fa/(q+1);
dfedlambda = dfedlambda  - q/(q+1) *v1 + ((1-lambda)*q+1)/(q+1) *dphi1dlam*(-log(1-phi1-alpha*phi1)/phi1^2 - (1+alpha)/phi1 - chi - cs/phi1^2 + cs/phi1^2 *(1+(alpha*phi1/cs)^2)^0.5);
dfedlambda = dfedlambda  + q/(q+1) *v2 + lambda*q/(q+1) *dphi2dlam*(-log(1-phi2-alpha*phi2)/phi2^2 - (1+alpha)/phi2 - chi - cs/phi2^2 + cs/phi2^2 *(1+(alpha*phi2/cs)^2)^0.5);

f(1) = dfedx;
f(2) = dfedlambda;

endfunction

//////////////////
// Main program //
//////////////////

//////////////////////
// Global variables //
//////////////////////
global q sig chi alpha cs

q = 4;
chi=0;
alpha=0;
cs=10^(-4);

y0 = [4.0; 0.35];

par0 = 0.08;
par9 = 0.36;
n=100;

//y0 = [6.21; 0.68]//
//
//par0 = 0.6;
//par9 = 0.26;
//n=200;

pst=(par9-par0)/n;
for k=1:n-1
  par=par0+pst*(k-1);// par=par0+pst*(k-1);
  sig=par;
    printf("\nsig=%f", sig);
  res(k,1)=par;
  [ymin,fval,info]=fsolve(y0,dfe);
  y0=ymin;                              //golden rule
  printf("\nforce=%f  lambda=%f   fval_1=%e  fval_2=%e\n", ymin(1), ymin(2), fval(1), fval(2));
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


s0 = '~/varvara/Documents/IMC/Analytical theory/Charged/max_ex/';
sout = s0 + 'lambda_vs_sig_q' + string(q) + '_v' + string(alpha) +'_max.dat';
u=file('open',sout,'unknown');
write(u, res, '(f12.6,32(e16.6))');
file('close',u);

// Plot
subplot(2,2,1); plot(res(:,1),res(:,2),'-r')
xlabel('sig')
ylabel('force')

subplot(2,2,2); plot(res(:,1),res(:,3),'-g')
xlabel('sig')
ylabel('lambda')

subplot(2,2,3); plot(res(:,1),res(:,4),'-b')
xlabel('sig')
ylabel('fval_1')

subplot(2,2,4); plot(res(:,1),res(:,5),'-m')
xlabel('sig')
ylabel('fval_2')



// Plot
//subplot(4,4,1); plot(res(:,1),res(:,2),'-r')
//xlabel('sig')
//ylabel('force_1')
////
//subplot(4,4,2); plot(res(:,1),res(:,3),'-r')
//xlabel('sig')
//ylabel('lambda')
////
//subplot(4,4,3); plot(res(:,1),res(:,4),'-r',res(:,1),res(:,5),'-g')
//xlabel('sig')
//ylabel('fval')
////
//// Plot
//subplot(4,4,4); plot(res(:,1),res(:,16),'-r')
//xlabel('sig')
//ylabel('force_conf.star')
////
//// Plot
//subplot(4,4,5); plot(res(:,1),res(:,17),'-r')
//xlabel('sig')
//ylabel('force_sol_1')
////
//// Plot
//subplot(4,4,6); plot(res(:,1),res(:,18),'-r')
//xlabel('sig')
//ylabel('force_sol_2')
////
//// Plot
//subplot(4,4,7); plot(res(:,1),res(:,19),'-r')
//xlabel('sig')
//ylabel('force_int_1')
////
//// Plot
//subplot(4,4,8); plot(res(:,1),res(:,20),'-r')
//xlabel('sig')
//ylabel('force_int_2')
////
//// Plot
//subplot(4,4,9); plot(res(:,1),res(:,21),'-r')
//xlabel('sig')
//ylabel('force_ch_1')
////
//// Plot
//subplot(4,4,10); plot(res(:,1),res(:,22),'-r')
//xlabel('sig')
//ylabel('force_ch_2')
////
//// Plot
//subplot(4,4,11); plot(res(:,1),res(:,23),'-r')
//xlabel('sig')
//ylabel('force_os.pr_1')
////
//// Plot
//subplot(4,4,12); plot(res(:,1),res(:,24),'-r')
//xlabel('sig')
//ylabel('force_os.pr_2')
//
