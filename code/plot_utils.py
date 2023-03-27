import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as feature
import numpy as np
import matplotlib.pyplot as plt
import shapefile as shp
#import seaborn as sns
import utils
import os
import pandas as pd
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
import cmasher as cmr

# Define the land feature
ccrs_land = feature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='black',
                                        facecolor='silver',
                                        alpha=1,
                                        linewidth=0.2)

# Define a the color palette for each sat
# colors = sns.color_palette("bright", len(utils.sats_dic))
# color_dic = dict(zip(list(utils.sats_dic.keys()), colors))


def generate_figure(coords, fig_size=20, dpi=200):

    """
    Generate a figure with a Robinson projection and a land feature on it at a given region.
    """

    plt.rcParams['font.family'] = 'monospace'

    ratio = 1.25
    fs = np.sqrt(fig_size)

    fig = plt.figure(figsize=(np.sqrt(ratio*fig_size),
                              np.sqrt(fig_size/ratio)),
                     dpi=dpi)

    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=0))
    ax.add_feature(ccrs_land, zorder=2)

    delta = 0.25
    coord_extent = [coords[0]-delta, coords[1] +
                    delta, coords[2]-delta, coords[3]+delta]
    ax.set_extent(coord_extent, crs=ccrs.PlateCarree())

    fig.tight_layout(h_pad=1)

    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      draw_labels=True,
                      linewidth=0.1*fs,
                      color='gray',
                      zorder=-2,
                      alpha=0.25,
                      linestyle='-')

    gl.top_labels = False
    gl.right_labels = False

    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()

    gl.xlabel_style = {'size': 1.35*fs}
    gl.ylabel_style = {'size': 1.35*fs}

    return fig, ax


def plot_date(ax, fig_size, day, month, year, x=1.115, y=0.98):

    # Plot the date on the right corner of the figure

    date = f'{int(day):02d}/{int(month):02d}/{int(year)}'
    ax.text(x, y,
            date,
            ha='left',
            va='top',
            #bbox=dict(facecolor='w', alpha=0.8, edgecolor='k'),
            transform=ax.transAxes,
            fontsize=1.45*np.sqrt(fig_size))


def plot_L4_field(datarray, ax, fig_size, cmap, vmin, vmax):

    # Plot the L4 field

    fs = np.sqrt(fig_size)

    nrt_l4 = datarray.plot.contourf(ax=ax,
                                    add_colorbar=False,
                                    cmap=cmap,
                                    transform=ccrs.PlateCarree(),
                                    linewidths=0.2*fs,
                                    vmax=vmax,
                                    vmin=vmin,
                                    levels=15,
                                    algorithm='mpl2005',
                                    zorder=-3)

    return nrt_l4


def plot_swot_orbit(ax, fig_size, lw=0.195, c='k', alpha=1):

    # Plot the SWOT orbit

    fs = np.sqrt(fig_size)

    file_nadir = '../data/swot_orbit/swot_calval_orbit_june2015-v2_nadir.shp'
    file_swath = '../data/swot_orbit/swot_calval_orbit_june2015-v2_swath.shp'

    sf_nadir = shp.Reader(file_nadir)
    sf_swath = shp.Reader(file_swath)

    lw = lw*fs
    for shape in sf_nadir.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        ax.plot(x, y, zorder=-2,
                transform=ccrs.PlateCarree(),
                color=c,
                alpha=alpha,
                linewidth=lw)

    for shape in sf_swath.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        ax.plot(x, y, zorder=-2,
                transform=ccrs.PlateCarree(),
                color=c,
                alpha=alpha,
                linewidth=lw)


def plot_path(ds,
              ax,
              fig_size,
              i,
              vmin,
              vmax,
              var='adt',
              c='red',
              ec='k',
              two_swaths=None,
              cmap='coolwarm',
              plot_geov=True):
    
    # Plot the ADT path and calculate the perpendicular geostrophic velocity

    fs = np.sqrt(fig_size)

    arr = np.array(ds['adt'])
    nans = np.isnan(arr)
    lws = np.full(arr.shape, 0.1*fs)
    lws[nans] = 0

    im = ds.plot.scatter(x='longitude',
                         y='latitude',
                         hue=var,
                         add_colorbar=False,
                         cmap=cmap,
                         transform=ccrs.PlateCarree(),
                         marker='o',
                         ec=ec,
                         lw=lws,
                         vmax=vmax,
                         vmin=vmin,
                         s=2.5*fs,
                         zorder=5)

    x = np.array(ds.longitude)
    y = np.array(ds.latitude)
    z = np.array(ds.geov)

    ax.plot(x, y,
            zorder=3,
            transform=ccrs.PlateCarree(),
            lw=0.1*fs,
            c=ec,
            linestyle='--')

    xs = x[0], x[-1]
    ys = y[0], y[-1]

    if y[-1] > y[0]:
        va = 'top'
        sgn = -1
    else:
        va = 'bottom'
        sgn = 1

    idx = 0

    if two_swaths == None:
        extra = ''
    else:
        extra = two_swaths

    ax.text(x[0], y[0]+0.015*fs*sgn,
            str(i)+extra,
            ha='center',
            va=va,
            transform=ccrs.PlateCarree(),
            fontsize=1.5*fs)

    def get_arrow(ds):

        lon = np.array(ds.longitude)
        lat = np.array(ds.latitude)

        dy = (lat[-1]-lat[0])
        dx = (lon[-1] - lon[0])
        norm = np.sqrt(dx**2+dy**2)
        dx = dx/norm
        dy = dy/norm

        return dx, dy

    dx, dy = get_arrow(ds)

    a = 0.1
    dx = a*dx
    dy = a*dy

    norm = np.sqrt(dx**2 + dy**2)

    angle = np.arctan(dy/dx)

    ax.arrow(x[-1]+0.02*fs*dx/norm,
             y[-1]+0.02*fs*dy/norm,
             dx,
             dy,
             width=0.00005*fs,
             head_width=0.0075*fs,
             head_length=0.01*fs,
             fc='black',
             ec='black',
             transform=ccrs.PlateCarree()
             )

    if plot_geov:
        q = ds.plot.quiver(x='longitude',
                           y='latitude',
                           u='geovx',
                           v='geovy',
                           pivot='tail',
                           width=0.0011*fs,
                           color=c,
                           ec='k',
                           headlength=2,
                           headaxislength=2,
                           headwidth=1.9,
                           lw=0.05*fs,
                           scale=5,
                           zorder=4,
                           add_guide=False,
                           transform=ccrs.PlateCarree())

        value = 0.1

        i = 0
        for value in [0.1, 0.25, 0.5]:

            qk = plt.quiverkey(q, 1.04, 0.31+i/30, value,
                               r'{:0.0f} cm s$^{{-1}}$'.format(value*100),
                               angle=45,
                               labelpos='W',
                               color='k',
                               fontproperties={'size': 1.2*fs},
                               coordinates='figure')
            i += 1

    return im


def plot_legend(sat_dic, ax, fig_size, swot=False, x=1.035, y=0.9):

    if swot == False:
        names = utils.sats_dic
        i = 1
    else:
        names = {'SWOT': 'SWOT'}
        i = 0

    labels = []
    times_ = []
    for sat in sat_dic:
        if not isinstance(sat_dic[sat], dict):
            t_i = pd.Timestamp(np.array(sat_dic[sat].time)[0])
            tt_i = f'{t_i.hour:02d}'+':'+f'{t_i.minute:02d}'
            labels.append(str(i)+'. ' + names[sat])
            times_.append(tt_i)

        else:
            times = []
            for swath in sat_dic[sat]:
                t_i = pd.Timestamp(np.array(sat_dic[sat][swath].time)[0])
                tt_i = f'{t_i.hour:02d}'+':'+f'{t_i.minute:02d}'
                times.append(tt_i)
            tt_i = '\n'.join(times)
            labels.append(str(i)+'. ' + names[sat]+'\n')
            times_.append(tt_i)

        i += 1

    label = '\n'.join(labels)
    time = '\n'.join(times_)

    ax.text(x, y,
            label,
            ha='left',
            va='top',
            #bbox=dict(facecolor='w', alpha=0.8, edgecolor='k'),
            transform=ax.transAxes,
            fontsize=1.15*np.sqrt(fig_size),
            family='monospace')

    ax.text(x + 0.22, y,
            time,
            ha='left',
            va='top',
            #bbox=dict(facecolor='w', alpha=0.8, edgecolor='k'),
            transform=ax.transAxes,
            fontsize=1.15*np.sqrt(fig_size),
            family='monospace')


def plot_along_track_sats(sat_dic,
                          ax,
                          fig_size,
                          vmin,
                          vmax,
                          cmap,
                          c=None,
                          i=1,
                          plot_geov=True):
    
    """
    Given a dictionary of satellite data, plot the along track data
    """

    for sat in sat_dic:

        if c == None:
            c = color_dic[sat]

        if not isinstance(sat_dic[sat], dict):
            a = 1

            im = plot_path(sat_dic[sat],
                           ax,
                           fig_size,
                           i,
                           vmin,
                           vmax,
                           c=c,
                           cmap=cmap,
                           plot_geov=plot_geov)

        else:
            j = 0
            subindex = ['\u2081', '\u2082']
            for swath in sat_dic[sat]:
                im = plot_path(sat_dic[sat][swath],
                               ax,
                               fig_size,
                               i,
                               vmin,
                               vmax,
                               c=c,
                               two_swaths=subindex[j],
                               cmap=cmap,
                               plot_geov=plot_geov)
                j += 1

        i += 1

    return im


def plot_cbar(im,
              fig, ax,
              fig_size,
              label='ADT (m)',
              location='left',
              extend='both'):

    fs = np.sqrt(fig_size)

    rot = 90
    scalar = 0.01
    anchor = (1.2, 0.95)
    if location == 'bottom':
        rot = 0
        scalar = 0.02
        anchor = (0.5, 1.2)

    cbar = fig.colorbar(im,
                        location=location,
                        extend=extend,
                        shrink=0.75,
                        aspect=30,
                        pad=scalar*fs,
                        anchor=anchor,
                        ax=ax)

    cbar.ax.xaxis.set_major_formatter('{:.2f}'.format)
    cbar.ax.yaxis.set_major_formatter('{:.2f}'.format)
    cbar.ax.tick_params(labelsize=1.5*fs)
    cbar.set_label(label, rotation=rot,
                   fontsize=1.5*fs,
                   # labelpad=lp*fig_size
                   )


def get_map(ds_nrt, swot, sat_dic, product, day, month, year, coords, dpi=250, save=True):

    """
    Plot all data (L4, L3 and SWOT-nadir)
    """

    vmax_l3, vmin_l3, vmax_l4, vmin_l4 = utils.get_min_max(
        ds_nrt, sat_dic, swot, product)

    if product == 'altimetry':
        title = 'ADT and along track \u22A5 geostrophic velocities'
        cbar_label = 'ADT (m)'
        cmap_l3 = cmr.get_sub_cmap('jet', 0, 1)
        cmap_l4 = cmr.get_sub_cmap('jet', 0, 1)
        arrow_color = 'fuchsia'

    elif product == 'sst':
        title = 'SST and along track ADT and \u22A5 geostrophic velocities'
        cbar_label = 'SST ($^\circ \!$C)'
        cmap_l3 = cmr.get_sub_cmap('jet', 0, 1)
        cmap_l4 = cmr.get_sub_cmap('coolwarm', 0.51, 1)
        arrow_color = 'c'

    elif product == 'chlorophyll':
        title = 'CHL and along track ADT and \u22A5 geostrophic velocities'
        cbar_label = 'CHL (mg m$^{-3}$)'
        cmap_l3 = cmr.get_sub_cmap('jet', 0, 1)
        cmap_l4 = cmr.get_sub_cmap('summer', 0, 0.8)
        arrow_color = 'r'
        vmin_l4 = 0

    fig_size = 20
    fs = np.sqrt(fig_size)

    fig, ax = generate_figure(coords, fig_size=fig_size, dpi=250)
    plot_date(ax, fig_size, day, month, year)

    if swot:
        im_l3 = plot_along_track_sats(
            swot, ax, fig_size, vmin_l3, vmax_l3, cmap=cmap_l3, c='w', i=0, plot_geov=False)
        plot_legend(swot, ax, fig_size, swot=True)

    if sat_dic:
        im_l3 = plot_along_track_sats(
            sat_dic, ax, fig_size, vmin_l3, vmax_l3, cmap=cmap_l3, c=arrow_color)
        plot_legend(sat_dic, ax, fig_size, y=0.8)

    if swot or sat_dic:
        plot_cbar(im_l3, fig, ax, fig_size, extend='both')

    im_l4 = plot_L4_field(ds_nrt,
                          ax,
                          fig_size,
                          cmap=cmap_l4,
                          vmin=vmin_l4,
                          vmax=vmax_l4)

    plot_swot_orbit(ax, fig_size, c='w')

    plot_cbar(im_l4,
              fig, ax, fig_size,
              label=cbar_label, location='bottom', extend='neither')

    ax.set_title(title,
                 fontsize=1.5*fs)

    if save == True:
        savepath = f'../figures/{product}/{day}_{month}_{year}.png'

        dir_path = os.path.dirname(os.path.realpath(savepath))+'/'
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        plt.savefig(savepath, dpi=dpi, bbox_inches='tight')

    return plt.show()
