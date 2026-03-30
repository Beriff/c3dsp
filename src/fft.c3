module c3dsp;

import std::math::complex;
import std::math;
import std::io;


fn Complex complex_exp(double x) {
  return (Complex){.r = math::cos(x), .c = math::sin(x)};
}

fn Complex[] fft(Complex[] x) {
  int n = x.len;
	Complex[] out = mem::new_array(Complex,n);
  
  for(int k = 0; k < n; k++) {
    for(int i = 0; i < n;i++) 
      out[k] += x[i] * complex_exp(-2.0 * math::PI * k*i/n);
  } 
  
	return out;
}

