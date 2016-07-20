q=7;
f=q+1;
npp=100;
M=(%pi*npp/2)/acos(sqrt((f-1)/f));
for k=1:200
         V(k)= -3*log(cos((%pi/2)*(k/200)))-((3*%pi^2/8)*((k/200)^2-(k/M)^2));
end
for k=1:200
    res(k,1)=k;
    res(k,2)=V(k);
end // for k
//
// Write data into file
s0 = '~/imc/StarBrush/Neutral/Poty/test/' ;
s = s0 + 'q'+string(q)+'_theory_z_V*' + '.dat';
u=file('open',s,'unknown');
write(u, res, '(1(f14.8), 32(e16.8))');
file('close',u);
