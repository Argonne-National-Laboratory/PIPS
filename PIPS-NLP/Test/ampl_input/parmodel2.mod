param s;    
set P:={1..4};     
var x{P};


minimize obj:  (x[1]+x[2])^2 +  (x[1]+x[2])*x[3] + (x[1]+x[2]+x[3])*x[4];

s.t. c0: x[1]*x[2]*x[4] == 10;
s.t. c1: x[1]*x[3] + x[2]^2 <= 5;
