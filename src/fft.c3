module c3dsp;

import std::math::complex;
import std::math;


<*
 fft function calculate the discrete Fourier transform of x or the inverse(normalized) transform if inverse is true.

 Note that algorithm work most efficient when the x.len is a power of 2.
*>
fn Complex[] fft(Complex[] x, bool inverse)
{
	if (is_power_of_2(x.len))
	{
		Complex[] res = stockham_fft(x, (inverse) ? -1 : 1);
		if (inverse)
		{
			for (int i = 0; i < res.len; i++)
			{
				res[i] *= (Complex) { .r = 1.0 / res.len, .c = 0.0 };
			}
		}

		return res;
	}

	return (inverse) ? bluesteini_ftt(x) : bluestein_ftt(x);
}

fn Complex complex_exp(double x) @private
{
	return (Complex) { .r = math::cos(x), .c = math::sin(x) };
}

macro bool is_power_of_2(int n) @private
{
	return !(n & (n - 1) != 0 || n < 1);
}

<*
 Naive FT approach based on definition.
*>
fn Complex[] naive_ft(Complex[] x, bool inverse) @private
{
	int n = x.len;
	Complex[] out = mem::new_array(Complex, n);

	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			out[k] += (Complex) { .r = (inverse) ? (1.0 / n) : (1.0), .c = 0.0 } * x[i]
				* complex_exp(-2.0 * math::PI * k * i / n);
		}
	}

	return out;
}

<*
 Stockham fft is implementation of discrete Fourier transform. It is very efficient, but limited to x.len = 2^n.
 Pass dir = 1 for forward unnormalized and -1 for inverse unnormalized.
*>
fn Complex[] stockham_fft(Complex[] x, int dir)
{
	int n = x.len;
	if (!is_power_of_2(n))
	{
		return x;
	}

	Complex[] out = @clone_slice(x);

	Complex[] tmp = mem::new_array(Complex, n);
	defer
	{
		free(tmp);
	}
	for (int stage = (int)std::math::log2((double)n); stage >= 1; stage--)
	{
		int nhalf = n >> stage; // nhalf = n / 2^stage
		int stride = 1 << (stage - 1); // stride = 2^(stage-1)

		for (int p = 0; p < nhalf; p++)
		{
			// Twiddle factor. It is constant for this p (applied only to the "odd" leg)
			Complex w = complex_exp(-2.0 * std::math::PI * p / (2.0 * nhalf) * (double)dir);

			for (int q = 0; q < stride; q++)
			{
				int idx0_i = q + stride * (2 * p + 0);
				int idx1_i = q + stride * (2 * p + 1);

				// Stockham autosort re-indexing for output (this is what eliminates the bit-reversal pass)
				int idx0_o = q + stride * (p + 0);
				int idx1_o = q + stride * (p + nhalf);

				Complex a = out[idx0_i];
				Complex b = w * out[idx1_i];

				// Radix-2 DIT butterfly
				tmp[idx0_o] = a + b;
				tmp[idx1_o] = a - b;
			}
		}

		free(out);
		out = @clone_slice(tmp);
	}

	return out;
}

<*
 Bluestein algorithm for prime-n sized FFTs, i.e. pad input to a size m = 2^p >= 2^(n-1).
*>
fn Complex[] bluestein_ftt(Complex[] x)
{
	int n = x.len;
	if (n == 0)
	{
		return x;
	}

	// Compute m = smallest power of 2 >= 2*n - 1 (padding for convolution)
	int nn_val = 2 * n - 1;
	double log_val = std::math::log2((double)nn_val);
	int m = 1 << (int)std::math::floor(log_val);
	if (m < nn_val)
	{
		m <<= 1;
	}

	// y = padded input + chirp pre-multiplication
	Complex[] y = mem::new_array(Complex, m);
	for (int i = 0; i < n; i++)
	{
		y[i] = x[i];
	}
	defer
	{
		free(y);
	}

	// w = chirp filter (symmetric)
	Complex[] w = mem::new_array(Complex, m);
	w[0] = (Complex) { .r = 1.0, .c = 0.0 };
	defer
	{
		free(w);
	}

	double a0 = std::math::PI / (double)n;
	for (int i = 1; i < n; i++)
	{
		Complex angle = complex_exp(a0 * (double)i * (double)i);
		w[i] = w[m - i] = angle;

		// y[i] *= exp(-j * theta)
		y[i] = y[i] * (Complex) { .r = angle.r, .c = -angle.c };
	}

	// Forward FFT (unnormalized)
	Complex[] y_fft = stockham_fft(y, 1);
	free(y);
	y = y_fft;

	// FFT of the chirp filter
	Complex[] w_fft = stockham_fft(w, 1);
	defer
	{
		free(w_fft);
	}

	// Pointwise multiply (circular convolution)
	for (int i = 0; i < m; i++)
	{
		y[i] *= w_fft[i];
	}

	// Inverse FFT (unnormalized)
	Complex[] y_ifft = stockham_fft(y, -1);
	free(y);
	y = y_ifft;

	// Final scaling with conjugate chirp / m
	for (int i = 0; i < n; i++)
	{
		Complex scaled_conj = (Complex) { .r = w[i].r * 1.0 / (double)m, .c = -w[i].c * 1.0 / (double)m };
		y[i] *= scaled_conj;
	}

	// Return only first n elements
	Complex[] res = @clone_slice(y[0..n - 1]);

	return res;
}


<*
 Inverse bluestein ftt.
*>
fn Complex[] bluesteini_ftt(Complex[] x)
{
	int n = x.len;
	if (n == 0)
	{
		return x;
	}

	// Conjugate input
	Complex[] y = mem::new_array(Complex, n);
	defer
	{
		free(y);
	}
	for (int i = 0; i < n; i++)
	{
		y[i] = (Complex) { .r = x[i].r, .c = -x[i].c };
	}

	// Bluestein forward on conjugated input
	Complex[] res = bluestein_ftt(y);

	// Conjugate output and scale by 1/n
	for (int i = 0; i < n; i++)
	{
		res[i] = (Complex) { .r = res[i].r * 1.0 / (double)n, .c = -res.[i].c * 1 / (double)n };
	}

	return res;
}
