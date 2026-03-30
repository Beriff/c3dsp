module c3dsp;
import std::io;
import std::math::complex;

fn int main(String[] args)
{
	io::printn(fft({{0,0},{1,0},{2,0}}));
	return 0;
}
