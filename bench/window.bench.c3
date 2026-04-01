module c3dsp::bench @benchmark;

import c3dsp;
import std::math;

fn void window_bench_rect() {
	Complex[] w = c3dsp::get_window(WindowType.RECTANGULAR, WindowMode.SYMMETRIC, 2048);

	free(w);
}

fn void window_bench_blackmanharris() {
	Complex[] w = c3dsp::get_window(WindowType.BLACKMANHARRIS, WindowMode.SYMMETRIC, 2048);

	free(w);
}