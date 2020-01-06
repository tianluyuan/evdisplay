# evdisplay
`chmod a+x evdisplay.py`

Full event
`./evdisplay.py -p InIcePulses_BaselineCorrected -d --tlim 13300 18580 --xlim 150 450 --ylim 200 500 --zlim -400 -200 -i $I3_GCD ~/Downloads/Monopodrec_Hydrangea_more_iter_newbaseline.i3 --view 30 108`

Early with cherenkov
`./evdisplay.py -p InIcePulses_BaselineCorrected -d --tlim 13300 13580 --xlim 150 450 --ylim 200 500 --zlim -400 -200 -i $I3_GCD ~/Downloads/Monopodrec_Hydrangea_more_iter_newbaseline.i3 --particle 511.076 349.199 -386.932 13253 1.2597 -0.1487 2 --scaling 0.1 --view 30 108 --cherenkov`
