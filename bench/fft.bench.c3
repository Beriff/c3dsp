module c3dsp::bench @benchmark;

import c3dsp;
import std::math;

fn void fft_bench_short() {
	Complex[] signal = { { 1, 0 }, { 2, 1 }, { 3, -1 }, { 4, 0 } };
	Complex[] res = c3dsp::fft(signal, false);
	free(res);
}

fn void fft_bench_long() {
	Complex[] signal = { 
		{ 1, 0 }, { 2, 1 }, { 3, -1 }, { 4, 0 },
		{ 1, 0 }, { 2, 1 }, { 3, -1 }, { 4, 0 },
		{ 1, 0 }, { 2, 1 }, { 3, -1 }, { 4, 0 },
		{ 1, 0 }, { 2, 1 }, { 3, -1 }, { 4, 0 },
		{ 1, 0 }, { 2, 1 }, { 3, -1 }, { 4, 0 },
	};

	Complex[] res = c3dsp::fft(signal, false);
	free(res);
}