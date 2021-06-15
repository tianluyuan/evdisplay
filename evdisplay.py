#!/usr/bin/env python3
import os
import argparse
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from collections import defaultdict, namedtuple
from I3Tray import I3Tray, I3Units
from icecube import icetray, dataio, dataclasses, phys_services


Particle = namedtuple('Particle', 'x y z t zen azi topo')
NIDX = 1.366 # n_g at 400 nm
C = 3e8 * I3Units.m/I3Units.s #m/s


def dir_to_mom(zen, azi):
    return np.asarray([-np.sin(zen)*np.cos(azi), -np.sin(zen)*np.sin(azi), -np.cos(zen)])


def range3d(x, p, step, stop):
    steps = np.arange(0,stop+step,step)
    steps[-1] = stop
    return np.asarray(p)[:,None]*steps[None,:]+np.asarray(x)[:,None]


def cend(particle, tstop):
    delta = (tstop-particle.t) * C
    p = dir_to_mom(particle.zen, particle.azi)
    x = np.asarray(particle[:3])
    return p*delta+x


def particle_nodes(particle, tstop, step):
    delta = (tstop-particle.t) * C
    p = dir_to_mom(particle.zen, particle.azi)
    x = np.asarray(particle[:3])
    if particle.topo == 0:
        yield range3d(x, p, step, min(20, delta))
    elif particle.topo == 1:
        yield range3d(x, p, step, delta)
    elif particle.topo == 2:
        sinfo = [_ for _ in particle]
        sinfo[-1] = 0
        shower_nodes = particle_nodes(Particle(*sinfo), tstop, step).next()
        yield shower_nodes
        _x, _y, _z = shower_nodes[:,-1]
        track = Particle(_x, _y, _z, min(particle.t+20./C, tstop),
                         particle.zen, particle.azi, 1)
        track_nodes = particle_nodes(track, tstop, step).next()
        yield track_nodes


def detector(i3_omgeo, all_pulses):
    oms = defaultdict(lambda: defaultdict(list))
    for om, omgeo in i3_omgeo:
        # if all_pulses.has_key(om):
        #     continue
        oms[omgeo.omtype][om.string].append(omgeo.position)
    return oms
        

def event(i3_omgeo, all_pulses, tlim):
    eve = defaultdict(list)
    for om in all_pulses.keys():
        om_pulses = all_pulses[om]
        fhit = np.inf
        qtot = 0
        for pulse in om_pulses:
            if pulse.time > tlim[1]:
                break
            if tlim[0] <= pulse.time:
                qtot += pulse.charge
                fhit = min(fhit, pulse.time)
        if qtot > 0:
            eve['q'].append(qtot)
            eve['t'].append(fhit)
            eve['x'].append(i3_omgeo[om].position.x)
            eve['y'].append(i3_omgeo[om].position.y)
            eve['z'].append(i3_omgeo[om].position.z)
    return eve


def draw_sphere(center, r, ax, **kwargs):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = r * np.outer(np.cos(u), np.sin(v)) + center[0]
    y = r * np.outer(np.sin(u), np.sin(v)) + center[1]
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
    ax.plot_surface(x, y, z, **kwargs)


def draw_cherenkov(nodes, stop, ax, **kwargs):
    for _ in nodes.T:
        r = np.linalg.norm(stop-_)/NIDX
        draw_sphere(_, r, ax, **kwargs)


def colorlines3d(x, y, z, ax, ncolors=5, cmapname='viridis_r', **kwargs):
    """Plot a line plot in which the lines change colors as the data is
    stepped through.

    *ncolors* specifies the number of different colors to use
    """
    cmap = plt.get_cmap(cmapname)
    norm = colors.Normalize(vmin=0, vmax=ncolors-1)
    for i in range(ncolors):
        chunksize = len(x)//ncolors
        low = i*chunksize
        # add 1 to keep lines connected
        high = min((i+1)*chunksize+1, len(x))
        ax.plot(x[low:high], y[low:high], z[low:high], color=cmap(norm(i)), **kwargs)

        
def draw_event(frame, draw_detector, draw_coord, draw_grid, pulse,
               xlim, ylim, zlim, tlim, particle, step,
               scaling, cmap, depthshade, cherenkov,
               view, llhout):
    """ plots the pulses of the individual frame
    """
    if not frame.Has(pulse):
        return

    run_id = frame['I3EventHeader'].run_id
    event_id = frame['I3EventHeader'].event_id
    mjd = np.round(frame['I3EventHeader'].start_time.mod_julian_day_double, 2)

    all_pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(
        frame, pulse)
    dom_cal = frame['I3Calibration'].dom_cal
    i3_omgeo = frame['I3Geometry'].omgeo

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if draw_grid:
        X, Y = np.meshgrid(np.arange(-600, 650, 50), np.arange(-600, 650, 50))
        ax.plot_wireframe(X, Y, np.full(Y.shape,-500), color='gray', linewidth=0.1)
    if draw_detector:
        oms = detector(i3_omgeo, all_pulses)
        for omtype in oms:
            for s in oms[omtype]:
                _ = np.asarray(oms[omtype][s]).T
                if omtype in [dataclasses.I3OMGeo.OMType.IceCube,
                              dataclasses.I3OMGeo.OMType.PDOM,
                              dataclasses.I3OMGeo.OMType.DEgg,
                              dataclasses.I3OMGeo.OMType.mDOM]:
                    ax.plot(*_, 'k-', linewidth=0.3)
                ax.scatter(*_, marker='.',
                           edgecolor='none',s=3, c='k', depthshade=depthshade)
        # for _x, _y in zip(str_x, str_y):
        #     ax.plot([_x]*2, [_y]*2, [om_z[0], om_z[-1]], c='grey')

    if None not in particle:
        # get the distance travelled by c from vertex to tlim[1]
        stop = cend(particle, tlim[1])

        for _ in particle_nodes(particle, tlim[1], step):
            pobj = ax.plot(*_)
            if cherenkov:
                draw_cherenkov(_, stop, ax,
                               color=pobj[-1].get_color(), alpha=0.1)
    if os.path.isfile(llhout):
        llh = pd.read_csv(
            llhout, delim_whitespace=True, header=None,
            names='l rlogl x y z zenith azimuth e t a b'.split(),
            error_bad_lines=False)
        llhsteps = llh.loc[llh['l'].str.isdigit()].apply(pd.to_numeric)
        colorlines3d(llhsteps['x'], llhsteps['y'], llhsteps['z'], ax)
        
    eve = event(i3_omgeo, all_pulses, tlim)
    ax.scatter(eve['x'], eve['y'], eve['z'], marker='.', edgecolor='none', s=np.asarray(eve['q'])/scaling,
               c=np.log(np.asarray(eve['t'])-min(eve['t'])), cmap=cmap, depthshade=depthshade)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    if not draw_coord:
        ax.axis("off")
    ax.view_init(*view)
    plt.show()


def main():
    """ plot steamshovel-like eventdisplay with matplotlib
    """
    parser = argparse.ArgumentParser(
        description='Input i3 file and output matplotlib evdisplay')
    parser.add_argument('-i', '--infile', nargs='+')
    parser.add_argument('-g', '--grid', default=False, action='store_true')
    parser.add_argument('-d', '--det', default=False, action='store_true')
    parser.add_argument('-c', '--coord', default=False, action='store_true')
    parser.add_argument('-N', '--nframes', type=int, default=None, help='number of frames to process')
    parser.add_argument('-p', '--pulse', type=str,
                        default='OfflinePulsesHLC',
                        help='specify the pulse type (default SplitInIcePulses)')
    parser.add_argument('--xlim', nargs=2, default=(None,None), type=float, help='(xlow, xup)')
    parser.add_argument('--ylim', nargs=2, default=(None,None), type=float, help='(ylow, yup)')
    parser.add_argument('--zlim', nargs=2, default=(-600, 500), type=float, help='(zlow, zup)')
    parser.add_argument('--tlim', nargs=2, default=(-np.inf, np.inf), type=float, help='(tlow, tup)')
    parser.add_argument('--particle', nargs=7, default=(None, None, None, None, None, None, None),
                        type=float, help='(x,y,z,t,zen,azi,topo) topo=0 cascade, 1 track, 2 hybrid')
    parser.add_argument('--step', default=10, type=float, help='nodes for particle, cherenkov bubbles placed here')
    parser.add_argument('-s', '--scaling', default=1,
                        type=float, help='factor to scale down qtot by for bubble size')
    parser.add_argument('--cmap', default='jet_r')
    parser.add_argument('--depthshade', default=False, action='store_true')
    parser.add_argument('--cherenkov', default=False, action='store_true', help='draw cherenkov sphere')
    parser.add_argument('--view', nargs=2, default=(None,None), type=float, help='Passed to ax.view_init for initial elev and azimuth angle')
    parser.add_argument('--llhout', default=None, help='Path to DF llh 10 output file to show vertex convergence')

    args = parser.parse_args()

    tray = I3Tray()
    tray.Add('I3Reader', Filenamelist=args.infile)
    tray.Add(draw_event,
             draw_grid=args.grid,
             draw_detector=args.det,
             draw_coord=args.coord,
             pulse=args.pulse,
             xlim=args.xlim,
             ylim=args.ylim,
             zlim=args.zlim,
             tlim=args.tlim,
             particle=Particle(*args.particle),
             scaling=args.scaling,
             cmap=args.cmap,
             depthshade=args.depthshade,
             cherenkov=args.cherenkov,
             step=args.step,
             view=args.view,
             llhout=args.llhout)

    if args.nframes is None:
        tray.Execute()
    else:
        tray.Execute(args.nframes)
    tray.Finish()

if __name__ == '__main__':
    main()
