module c3dsp;

import std::math::complex;
import std::io;

enum WindowType : int
{
	RECTANGULAR,
	HANN,
	HAMMING,
	BLACKMAN,
	BLACKMANHARRIS
	// TODO:
	// KAISER
	// FLATTOP
	// ...
}

enum WindowMode : int
{
	PERIODIC,
	SYMMETRIC
}

<*
	Return the coefficients of window function of length n. 
	If mode == PERIODIC, the periodic form of the window is returned(This is equivalent to the window of length n+1 with the last coefficient removed). 
*>
fn Complex[] get_window(WindowType type, WindowMode mode, int n)
{
	Complex[] w;

	double denom = 1.0;
	switch (mode)
	{
		case WindowMode.PERIODIC:
			denom = (double)n;
		case WindowMode.SYMMETRIC:
			denom = (double)(n - 1);
	}

	switch (type)
	{
		case WindowType.RECTANGULAR:
			w = rect_window(n);
		case WindowType.HANN:
			w = hann_window(n, denom);
		case WindowType.HAMMING:
			w = hamming_window(n, denom);
		case WindowType.BLACKMAN:
			w = blackman_window(n, denom);
		case WindowType.BLACKMANHARRIS:
			w = blackmanharris_window(n, denom);
	}

	return w;
}

<*
	Calculate the coefficients of window function and multiply by the given signal(x[i] *= w[i]). 
	If mode == PERIODIC, the signal multiply by the periodic form of the window. 
*>
fn void apply_window(WindowType type, WindowMode mode, Complex[]* x, int n)
{
	Complex[] w = get_window(type, mode, n);
	defer free(w);

	foreach (idx, &val : x) *val *= w[idx];
}

macro Complex[] rect_window(int n) @private
{
	Complex[] w = mem::new_array(Complex, n);
	foreach (idx, &val : w) *val = (idx <= (n - 1)) ? { 1.0, 0 } : { 0.0, 0 };

	return w;
}

macro Complex[] hann_window(int n, double denom) @private
{
	Complex[] w = mem::new_array(Complex, n);
	foreach (idx, &val : w) if (idx <= (n
		- 1)) *val = (Complex) { .r = 0.5 * (1 - math::cos(2 * math::PI * idx / denom)), .c = 0.0 };

	return w;
}

macro Complex[] hamming_window(int n, double denom) @private
{
	Complex[] w = mem::new_array(Complex, n);
	foreach (idx, &val : w) if (idx <= (n
		- 1)) *val = (Complex) { .r = 0.54 - 0.46 * math::cos(2 * math::PI * idx / denom), .c = 0.0 };

	return w;
}

macro Complex[] blackman_window(int n, double denom) @private
{
	Complex[] w = mem::new_array(Complex, n);
	foreach (idx, &val : w)
	{
		if (idx <= (n
			- 1)) *val = (Complex) {
				.r = 0.42 - 0.5 * math::cos(2 * math::PI * idx / denom) + 0.08 * math::cos(4 * math::PI * idx / denom),
				.c = 0.0
			};
	}

	return w;
}

macro Complex[] blackmanharris_window(int n, double denom) @private
{
	Complex[] w = mem::new_array(Complex, n);
	foreach (idx, &val : w)
	{
		if (idx <= (n
			- 1)) *val = (Complex) {
				.r = 0.35875 - 0.48829 * math::cos(2 * math::PI * idx / denom)
					+ 0.14128 * math::cos(4 * math::PI * idx / denom)
					- 0.01168 * math::cos(6 * math::PI * idx / denom),
				.c = 0.0
			};
	}

	return w;
}
