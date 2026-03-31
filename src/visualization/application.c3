module c3dsp::visualize;

import c3dsp;
import raylib55::rl;

struct VisualizationApp {
	int screenWidth;
	int screenHeight;
}

fn void VisualizationApp.new(self, int width, int height) {
	rl::set_trace_log_level(RLTraceLogLevel.WARNING);

	rl::init_window(width, height, "C3DSP Visualizer");
	rl::set_target_fps(60);
}

fn void VisualizationApp.beginUpdating(self) {
	while(!rl::window_should_close()) {
		rl::clear_background(rl::BLACK);

		rl::begin_drawing();
			rl::draw_fps(0,0);
		rl::end_drawing();
	}
	rl::close_window();
}