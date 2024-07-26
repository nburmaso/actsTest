#./build/actsTest mpd mpd 100000
#./build/actsTest mpd_noframe mpd_noframe 100000
#./build/actsTest mpd_notpc mpd_notpc 100000


if false; then
root -b -q 'analyse_seed_efficiency.C("mpd", 1.6)'
root -b -q 'analyse_seed_efficiency.C("mpd", 1.9)'
root -b -q 'analyse_seed_efficiency.C("mpd", 2.2)'
root -b -q 'draw_seed_efficiency.C("pi", "mpd","mpd","mpd")'
root -b -q 'draw_seed_efficiency.C("pr", "mpd","mpd","mpd")'
root -b -q 'analyse_tracking_efficiency.C("mpd", 1.6, 0.03)'
root -b -q 'analyse_tracking_efficiency.C("mpd", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("mpd", 2.2)'
root 'draw_tracking_efficiency.C("pi","mpd","mpd","mpd")'
root 'draw_tracking_efficiency.C("pr","mpd","mpd","mpd")'
fi


if false; then
root -b -q 'analyse_seed_efficiency.C("mpd_noframe", 1.6)'
root -b -q 'analyse_seed_efficiency.C("mpd_noframe", 1.9)'
root -b -q 'analyse_seed_efficiency.C("mpd_noframe", 2.2)'
root -b -q 'draw_seed_efficiency.C("pi", "mpd_noframe","mpd_noframe","mpd_noframe")'
root -b -q 'draw_seed_efficiency.C("pr", "mpd_noframe","mpd_noframe","mpd_noframe")'
root -b -q 'analyse_tracking_efficiency.C("mpd_noframe", 1.6, 0.03)'
root -b -q 'analyse_tracking_efficiency.C("mpd_noframe", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("mpd_noframe", 2.2)'
root 'draw_tracking_efficiency.C("pi","mpd_noframe","mpd_noframe","mpd_noframe")'
root 'draw_tracking_efficiency.C("pr","mpd_noframe","mpd_noframe","mpd_noframe")'
fi



if false; then
root -b -q 'analyse_seed_efficiency.C("mpd_notpc", 1.6)'
root -b -q 'analyse_seed_efficiency.C("mpd_notpc", 1.9)'
root -b -q 'analyse_seed_efficiency.C("mpd_notpc", 2.2)'
root -b -q 'analyse_tracking_efficiency.C("mpd_notpc", 1.6, 0.03)'
root -b -q 'analyse_tracking_efficiency.C("mpd_notpc", 1.9)'
root -b -q 'analyse_tracking_efficiency.C("mpd_notpc", 2.2)'
root -b -q 'draw_seed_efficiency.C("pi", "mpd_notpc","mpd_notpc","mpd_notpc")'
root -b -q 'draw_seed_efficiency.C("pr", "mpd_notpc","mpd_notpc","mpd_notpc")'
fi
root 'draw_tracking_efficiency.C("pi","mpd_notpc","mpd_notpc","mpd_notpc")'
root 'draw_tracking_efficiency.C("pr","mpd_notpc","mpd_notpc","mpd_notpc")'

