#ifndef tracker_config
#define tracker_config

std::vector<double> positions{{2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000}}; // z-positions of sensitive layers
double thickness{1_um};
double rMin{357_mm};
double rMax{1300_mm};
double bz = 0.5_T;

#endif