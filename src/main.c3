module c3dsp;

import c3dsp::visualize;
import std::io;
import std::math::complex;

fn int main(String[] args)
{
	Complex[] signal = { { 1, 0 }, { 2, 1 }, { 3, -1 }, { 4, 0 } };
	Complex[] expected = { { 10, 0 }, { -1, 3 }, { -2, -2 }, { -3, -1 } };
	Complex[] res = fft(signal, false);

	io::print("Result: ");
	io::print(res);
	io::print("\nExpected: ");
	io::printn(expected);

	free(res);

	visualize::VisualizationApp app = {};
	app.new(800,600);
	app.beginUpdating();

	return 0;
}
