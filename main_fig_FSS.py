import sys; sys.path.append("../../")
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


def save_fig(fig_name='test', fig_format='png'):
    if fig_format == 'png':
        plt.savefig(fig_name+'.png', format='png', dpi=600)
    elif fig_format == 'pdf':
        plt.savefig(fig_name+'.pdf', format='pdf', dpi=600)
    elif fig_format == 'svg':
        plt.savefig(fig_name+'.svg', format='svg', dpi=600)
    else:
        plt.show()


def set_style_PR(fig_width=3.25, aspect_ratio = 0.6):
    fig_height = aspect_ratio*fig_width
    params = {
        'figure.figsize': (fig_width,fig_height),
        'legend.fontsize': 6,
        'legend.frameon': False,
        'axes.labelsize': 7,
        'axes.linewidth': 1.,
        'axes.linewidth': 0.8,
        'lines.linewidth': 0.75,
        'lines.markersize': 1.5,
        'xtick.labelsize' :7,
        'ytick.labelsize': 7,
        'mathtext.fontset': 'stixsans',
        'mathtext.rm': 'serif',
        'mathtext.bf': 'serif:bold',
        'mathtext.it': 'serif:italic',
        'mathtext.sf': 'sans\\-serif',
        'font.size':  7,
        'font.family': 'serif',
        'font.serif': "Helvetica",
    }
    mpl.rcParams.update(params)


def figure_02(o_name='fig', o_format='png'):

    # -- FIGURE LAYOUT --------------------------------------------------------

    set_style_PR(3.4, 0.45)
    fig = plt.figure()
    plt.subplots_adjust(left = 0.1, bottom = 0.18, right = 0.98, top = 0.97)
    gs00 = GridSpec(nrows = 1, ncols = 1)

    gsA = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs00[0,0], wspace=0.35, hspace=0.1)
    ax1 = fig.add_subplot(gsA[0, 0])
    ax2 = fig.add_subplot(gsA[0, 1])


    # -- CONVENTIENT FUNCTIONS AND OTHER DEFINITIONS --------------------------

    def subfig_label(ax, label):
        pos = ax.get_position()
        fig.text(
            pos.x0,
            pos.y1,
            label,
            fontsize=7,
            color="white",
            backgroundcolor="k",
            bbox=dict(facecolor="k", edgecolor="none", boxstyle="square,pad=0.1"),
            verticalalignment="top",
            horizontalalignment="left",
        )

    def fetch_data(f_name):
        dat = np.loadtxt(f_name)
        return dat[:,0], dat[:,1], dat[:,2]

    L_vals = [16,32,64,128,256,512]
    markers = ['o', 's', '>', '<', 'D', 'h']
    cols = plt.cm.cividis_r(np.linspace(0.1,1.,len(L_vals)))
    #cols = plt.cm.Blues(np.linspace(0.5,1.,len(L_vals)))

    xc=0.59269
    a=0.748
    b=0.105
    _sx = lambda x, L: (x-xc)*L**a
    _sy = lambda y, L: y*L**b

    # -- DATA CURVES FOR SUBFIGURES (A,B) -------------------------------------

    L = 16
    for idx, L in enumerate(L_vals):
        x, y, y_err = fetch_data('./ORDER_PARAMETER/orderParam_L%d.dat'%(L))
        ax1.plot(x, y, color=cols[idx], marker=markers[idx], label=r'$L=%d$'%(L))
        ax2.plot(_sx(x,L), _sy(y,L), color=cols[idx], marker=markers[idx], label=r'$L=%d$'%(L))



    # -- AXES DETAILS SUBFIGURE (A) -------------------------------------------

    x_lim = (0.4,0.7)
    x_ticks = (0.4,0.5,0.6,0.7)
    ax1.set_xlim(x_lim)
    ax1.set_xticks(x_ticks)
    ax1.tick_params(axis="x", length=2.5, pad=1, top=False)
    ax1.set_xlabel(r"Control parameter $x$",labelpad=1)

    y_lim = (0,0.7)
    y_ticks = (0,0.2,0.4,0.6)
    ax1.set_ylim(y_lim)
    ax1.set_yticks(y_ticks)
    ax1.tick_params(axis="y", length=2.0, pad=1)
    ax1.set_ylabel(r"Order parameter $y(x,L)$")

    legend = ax1.legend(
        ncol=1,
        handlelength=1.,
        borderpad=0.1,
        handletextpad=0.25,
        columnspacing=1.0,
        labelspacing=0.2,
        fontsize=5.0,
        title_fontsize=5.0,
        labelcolor="k",
        loc="center left",
    )

    subfig_label(ax1, '(a)')


    # -- AXES DETAILS SUBFIGURE (B) -------------------------------------------

    pos = ax2.get_position()
    fig.text(
        pos.x0+0.6*(pos.x1-pos.x0),
        pos.y0+0.3*(pos.y1-pos.y0),
        "$x_c=%6.5lf$\n$a=%4.3lf$\n$b=%4.3lf$"%(xc,a,b),
        fontsize=6,
        color="k",
        verticalalignment="bottom",
        horizontalalignment="left",
    )

    x_lim = (-3.2,3.2)
    x_ticks = (-3,-2,-1,0,1,2,3)
    ax2.set_xlim(x_lim)
    ax2.set_xticks(x_ticks)
    ax2.tick_params(axis="x", length=2.5, pad=1, top=False)
    ax2.set_xlabel(r"Scaling variable $(x-x_c)L^a$",labelpad=1)

    y_lim = (0,1.1)
    y_ticks = (0,0.2,0.4,0.6,0.8,1.0)
    ax2.set_ylim(y_lim)
    ax2.set_yticks(y_ticks)
    ax2.tick_params(axis="y", length=2.0, pad=1)
    ax2.set_ylabel(r"Scaled order par. $y(x,L) L^b$")

    legend = ax2.legend(
        ncol=1,
        handlelength=1.,
        borderpad=0.1,
        handletextpad=0.25,
        columnspacing=1.0,
        labelspacing=0.2,
        fontsize=5.0,
        title_fontsize=5.0,
        labelcolor="k",
        loc="center left",
    )

    subfig_label(ax2, '(b)')


    # -- SAVE FIGURE ---------------------------------------------------------- 
    save_fig(fig_name=o_name, fig_format=o_format )


def main():
    figure_02(o_name='fig_FSS')


main()
