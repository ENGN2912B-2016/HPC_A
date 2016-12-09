splot sprintf("timestep%i.dat",i) u 1:2:3:4 with points palette notitle 
i=i+1000
if(i < n) reread