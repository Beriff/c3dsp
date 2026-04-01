import c3dsp;

import std::math::complex;
import std::math;

import utils;

fn void test_fft_empty() @test
{
	Complex[] x = mem::new_array(Complex, 0);
	Complex[] res = c3dsp::fft(x, false);
	assert (res.len == 0);
	free(res);
}

fn void test_fft_all_zeros() @test
{
	int n = 32;
	Complex[] x = mem::new_array(Complex, n);

	Complex[] res = c3dsp::fft(x, false);
	for (int i = 0; i < n; i++)
	{
		assert (utils::approx_equal(res[i], { 0.0, 0.0 }, 1e-12));
	}

	free(x);
	free(res);
}

fn void test_fft_single_element() @test
{
	Complex[] x = { { 1.0, 0.0 } };
	Complex[] res = c3dsp::fft(x, false); // forward
	assert (utils::approx_equal(res[0], { 1.0, 0.0 }));
	free(res);

	Complex[] res_inv = c3dsp::fft(x, true); // inverse
	assert (utils::approx_equal(res_inv[0], { 1.0, 0.0 }));
	free(res_inv);
}

fn void test_fft_power_of_two() @test
{
	Complex[] x = { { 1.0, 0.0 }, { 0.0, 0.0 }, { -1.0, 0.0 }, { 0.0, 0.0 } };

	Complex[] expected = { { 0.0, 0.0 }, { 2.0, 0.0 }, { 0.0, 0.0 }, { 2.0, 0.0 } };

	Complex[] res = c3dsp::fft(x, false);
	assert (utils::arrays_approx_equal(res, expected, 1e-10));
	free(res);
}

fn void test_ifft_recover() @test
{
	Complex[] x = { { 1.0, 0.0 }, { 2.0, 1.0 }, { 3.0, -1.0 }, { 4.0, 0.0 } };

	Complex[] freq = c3dsp::fft(x, false); // forward
	Complex[] recovered = c3dsp::fft(freq, true); // inverse

	assert (utils::arrays_approx_equal(recovered, x, 1e-9));
	free(freq);
	free(recovered);
}

fn void test_bluestein_prime() @test
{
	Complex[] x = { { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 } };

	Complex[] expected = { { 5.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 } };

	Complex[] res = c3dsp::fft(x, false);
	assert (utils::arrays_approx_equal(res, expected, 1e-9));
	free(res);
}

fn void test_fft_cosine() @test
{
	int n = 16;
	Complex[] x = mem::new_array(Complex, n);

	for (int i = 0; i < n; i++)
	{
		x[i] = (Complex) { .r = std::math::cos(2.0 * std::math::PI * 3.0 * (double)i / (double)n), .c = 0.0 };
	}

	Complex[] res = c3dsp::fft(x, false);

	assert (utils::complex_magnitude(res[3]) > 7.9); // ~ N/2
	assert (utils::complex_magnitude(res[13]) > 7.9);

	free(x);
	free(res);
}
