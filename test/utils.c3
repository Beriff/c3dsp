import c3dsp;

import std::math::complex;
import std::math;

macro bool approx_equal(Complex a, Complex b, double eps = 1e-9)
{
	return complex_magnitude(a - b) <= eps;
}

macro bool arrays_approx_equal(Complex[] a, Complex[] b, double eps = 1e-9)
{
	if (a.len != b.len) return false;
	for (int i = 0; i < a.len; i++)
	{
		if (!approx_equal(a[i], b[i], eps)) return false;
	}
	return true;
}

macro complex_magnitude(Complex a)
{
	return std::math::sqrt(a.r * a.r + a.c * a.c);
}