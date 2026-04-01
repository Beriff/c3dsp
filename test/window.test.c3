import c3dsp;

import std::math::complex;
import std::math;
import std::io;
import utils;

fn void test_window_rect_coefficient() @test
{
	Complex[] w = c3dsp::get_window(WindowType.RECTANGULAR, WindowMode.SYMMETRIC, 4);
	Complex[] expected = { { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 } };

	assert ((utils::arrays_approx_equal(w, expected, 1e-5)));

	free(w);
}

fn void test_window_hann_coefficient() @test
{
	// Symmetric
	Complex[] w_sym = c3dsp::get_window(WindowType.HANN, WindowMode.SYMMETRIC, 4);
	Complex[] expected_sym = { { 0.0, 0.0 }, { 0.75, 0.0 }, { 0.75, 0.0 }, { 0.0, 0.0 } };

	assert ((utils::arrays_approx_equal(w_sym, expected_sym, 1e-5)));

	// Periodic
	Complex[] w_per = c3dsp::get_window(WindowType.HANN, WindowMode.PERIODIC, 4);
	Complex[] expected_per = { { 0.0, 0.0 }, { 0.5, 0.0 }, { 1.0, 0.0 }, { 0.5, 0.0 } };

	assert ((utils::arrays_approx_equal(w_per, expected_per, 1e-5)));

	free(w_sym);
	free(w_per);
}

fn void test_window_hamming_coefficient() @test
{
	// Symmetric
	Complex[] w_sym = c3dsp::get_window(WindowType.HAMMING, WindowMode.SYMMETRIC, 4);
	Complex[] expected_sym = { { 0.08, 0.0 }, { 0.77, 0.0 }, { 0.77, 0.0 }, { 0.08, 0.0 } };

	assert ((utils::arrays_approx_equal(w_sym, expected_sym, 1e-5)));

	// Periodic
	Complex[] w_per = c3dsp::get_window(WindowType.HAMMING, WindowMode.PERIODIC, 4);
	Complex[] expected_per = { { 0.08, 0.0 }, { 0.54, 0.0 }, { 1.0, 0.0 }, { 0.54, 0.0 } };

	assert ((utils::arrays_approx_equal(w_per, expected_per, 1e-5)));

	free(w_sym);
	free(w_per);
}

fn void test_window_blackman_coefficient() @test
{
	// Symmetric
	Complex[] w_sym = c3dsp::get_window(WindowType.BLACKMAN, WindowMode.SYMMETRIC, 4);
	Complex[] expected_sym = { { -1.3878e-17, 0.0 }, { 0.63, 0.0 }, { 0.63, 0.0 }, { -1.3878e-17, 0.0 } };

	assert ((utils::arrays_approx_equal(w_sym, expected_sym, 1e-10)));

	// Periodic
	Complex[] w_per = c3dsp::get_window(WindowType.BLACKMAN, WindowMode.PERIODIC, 4);
	Complex[] expected_per = { { -1.3878e-17, 0.0 }, { 0.34, 0.0 }, { 1.0, 0.0 }, { 0.34, 0.0 } };

	assert ((utils::arrays_approx_equal(w_per, expected_per, 1e-10)));

	free(w_sym);
	free(w_per);
}

fn void test_window_blackmanharris_coefficient() @test
{
	// Symmetric
	Complex[] w_sym = c3dsp::get_window(WindowType.BLACKMANHARRIS, WindowMode.SYMMETRIC, 4);
	Complex[] expected_sym = { { 6.0e-05, 0.0 }, { 5.2057e-01, 0.0 }, { 5.2058e-01, 0.0 }, { 6.0e-05, 0.0 } };

	assert ((utils::arrays_approx_equal(w_sym, expected_sym, 1e-5)));

	// Periodic
	Complex[] w_per = c3dsp::get_window(WindowType.BLACKMANHARRIS, WindowMode.PERIODIC, 4);
	Complex[] expected_per = { { 6.0e-05, 0.0 }, { 2.1747e-01, 0.0 }, { 1.0, 0.0 }, { 2.1747e-01, 0.0 } };

	assert ((utils::arrays_approx_equal(w_per, expected_per, 1e-5)));

	free(w_sym);
	free(w_per);
}


fn void test_window_single_signal() @test
{
	Complex[] signal = { { 1.0, 0.0 }, { 1.0, 0.0 } };
	Complex[] w = c3dsp::get_window(WindowType.BLACKMANHARRIS, WindowMode.SYMMETRIC, 4);
	foreach (idx, &val : signal) *val *= w[idx];

	assert (signal[0] == (Complex) { .r = w[0].r, .c = 0.0 } && signal[1] == (Complex) { .r = w[1].r, .c = 0.0 });

	free(w);
}
