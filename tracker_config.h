#ifndef tracker_config
#define tracker_config
//std::vector<double> positions{{2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000}}; // z-positions of sensitive layers
std::vector<double> positions{{2100, 2325, 2550, 2775, 3000}}; // z-positions of sensitive layers
//std::vector<double> positions{{2100, 2550, 3000}}; // z-positions of sensitive layers
double radLenSilicon = 9.370_cm;
double radLenAluminium = 8.897_cm;

// double thickness{0.00112*radLenSilicon};
double thickness{0.1_mm};
double rMin{357_mm};
double rMax{1300_mm};
double bz = 0.5_T;

#endif