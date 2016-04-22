param s;    
set P:={1..3};     
var x{P};
param d{1..2}; 


minimize obj:  x[1]*x[1] +  x[2]*x[2] + x[3]*x[3];

s.t. cap: x[2] + x[3] ==  d[s];
s.t. bal: x[1] - x[2] == 0;
