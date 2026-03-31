module c3dsp;

import std::math::complex;
import std::math;
import std::io;


fn Complex complex_exp(double x) @private
{
	return (Complex) { .r = math::cos(x), .c = math::sin(x) };
}

macro bool is_power_of_2(int n) @private
{
	return !(n & (n - 1) != 0 || n < 1);
}

fn Complex[] fft(Complex[] x)
{
    if(is_power_of_2(x.len)) {
        return stockham_fft(x,1);
    }

    return naive_ft(x);
}

<*
    Naive FT approach based on definition
*>
fn Complex[] naive_ft(Complex[] x) @private
{
	int n = x.len;
	Complex[] out = mem::new_array(Complex, n);

	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			out[k] += x[i] * complex_exp(-2.0 * math::PI * k * i / n);
		}
	}

	return out;
}

<*
 Stockham fft is implementation of discrete Fourier transform. It is very efficient, but limited to x.len = 2^n
 Pass dir = 1 for forward unnormalized and -1 for inverse unnormalized
*>
fn Complex[] stockham_fft(Complex[] x, int dir)
{
	int n = x.len;
	if (!is_power_of_2(n))
	{
		return x;
	}

	Complex[] out = mem::new_array(Complex, n);
	out = @clone_slice(x);

	Complex[] tmp = mem::new_array(Complex, n);
    defer {free(tmp);};
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

		out = @clone_slice(tmp);
	}

	return out;
}

