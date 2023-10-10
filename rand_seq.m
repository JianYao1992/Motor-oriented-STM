function f = rand_seq(len,n,m) % len: total FCSP number; n: total trial number; m: maximal FCSP number in each trial
f = int16(m*rand(1,n));   
s = len - (sum(f) - f(n));     
count = 0;
index = 0;
r = 0;
while s - m > 1e-18 || s < 1e-18
  while s < 1e-18
    index = int16(rand*(n-1) + 1);
    t = f(index);
    r = rand;
    f(index) = int16(t*r);       
    s = s + t - int16(t*r);
    count = count + 1;
  end
  while s - m > 1e-18
    index = int16(rand*(n-1) + 1);
    t = f(index);         
    r = rand;
    f(index) = t + int16((m-t)*r);
    s = s - int16((m-t)*r);  
    count = count + 1;
  end
end
f(n) = s;                
end