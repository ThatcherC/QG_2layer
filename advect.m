function qdot=advect(q,u,v,k,l)
  qdot=1i*k.*fft2(u.*q)+1i*l.*fft2(v.*q);
