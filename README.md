# evdisplay
`chmod a+x evdisplay.py`

Full event
`./evdisplay.py -p InIcePulses_BaselineCorrected -d --tlim 13300 18580 --xlim 150 450 --ylim 200 500 --zlim -400 -200 -i $I3_GCD ~/Downloads/Monopodrec_Hydrangea_more_iter_newbaseline.i3 --view 30 88`

Early with cherenkov
`./evdisplay.py -p InIcePulses_BaselineCorrected -d --tlim 13300 13580 --xlim 150 450 --ylim 200 500 --zlim -400 -200 -i $I3_GCD ~/Downloads/Monopodrec_Hydrangea_more_iter_newbaseline.i3 --particle 509.55 347.28 -387.17 13251.66 1.2814 -0.0473 2 --scaling 0.1 --view 30 88 --cherenkov`
