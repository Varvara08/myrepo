lambda=0.16666666667; // 1/6
np=100; // number of units in the brush
q=2;
f=q+1;
h2n=0.25;
M=(%pi*np/2)/acos(sqrt((f-1)/f));
//
// Считаем потенциал по формуле
for k=1:2*np+5
    if k < h2n*2*np then
         V(k)= -3*log(cos(h2n*%pi/2))-((3*%pi^2/8)*((h2n)^2-(h2n*2*np/M)^2))-(-3*log(cos((%pi/2)*(k/200)))-((3*%pi^2/8)*((k/200)^2-(k/M)^2)));
    else
        V(k)=0;
    end //for if
    w(k)=exp(-V(k)); // exp(-potential)
end
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
