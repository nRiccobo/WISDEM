import os
import matplotlib.pyplot as plt
import matplotlib.text as txt
try:
    import ruamel.yaml as ry
except:
    try:
        import ruamel_yaml as ry
    except:
        raise ImportError('No YAML package found')


# ---------------------
def load_yaml(fname_input):
    reader = ry.YAML(typ="safe", pure=True)
    with open(fname_input, "r", encoding="utf-8") as f:
        input_yaml = reader.load(f)
    return input_yaml


# ---------------------
def write_yaml(instance, foutput):
    # Write yaml with updated values
    yaml = ry.YAML()
    yaml.default_flow_style = None
    yaml.width = float("inf")
    yaml.indent(mapping=4, sequence=6, offset=3)
    yaml.allow_unicode = False
    with open(foutput, "w", encoding="utf-8") as f:
        yaml.dump(instance, f)

# ---------------------
def write_optimized_yaml(forig, yaml_opt, ftmp):
    # Read original yaml
    with open(forig, 'r') as fr:
        orig = fr.readlines()

    # Write optimized files
    fw = open(ftmp, 'w')
    flag_tow = False
    flag_mono = False
    flag_diam = False
    flag_thick = False
    for line in orig:
        if not flag_tow and line.find('tower:') >= 0:
            flag_tow = True
        if not flag_mono and line.find('monopile:') >= 0:
            flag_mono = True
        if not flag_diam and (flag_tow or flag_mono) and line.find('outer_diameter:') >= 0:
            flag_diam = True
        if not flag_thick and (flag_tow or flag_mono) and line.find('thickness:') >= 0:
            flag_thick = True

        if flag_tow and flag_diam and line.find('values:') >= 0:
            flag_diam = False
            val = yaml_opt['components']['tower']['outer_shape_bem']['outer_diameter']['values']
            tok = line.split(':')
            fw.write(f'{tok[0]}: {val}\n')
        elif flag_tow and flag_thick and line.find('values:') >= 0:
            flag_thick = False
            flag_tow = False
            val = yaml_opt['components']['tower']['internal_structure_2d_fem']['layers'][0]['thickness']['values']
            tok = line.split(':')
            fw.write(f'{tok[0]}: {val}\n')
        elif flag_mono and flag_diam and line.find('values:') >= 0:
            flag_diam = False
            val = yaml_opt['components']['monopile']['outer_shape_bem']['outer_diameter']['values']
            tok = line.split(':')
            fw.write(f'{tok[0]}: {val}\n')
        elif flag_mono and flag_thick and line.find('values:') >= 0:
            flag_thick = False
            flag_mono = False
            val = yaml_opt['components']['monopile']['internal_structure_2d_fem']['layers'][0]['thickness']['values']
            tok = line.split(':')
            fw.write(f'{tok[0]}: {val}\n')
        else:
            fw.write(line)

    fw.close()
# ---------------------

def init_fig_axis():
    fig = plt.figure(figsize=(12,9))
    ax  = fig.add_subplot(111)
    return fig, ax
# ---------------------

def save(fig, froot, mode='png'):
    assert mode in ['png','eps','pdf','all']
    fileName, fileExtension = os.path.splitext(froot)
    padding = 0.1
    dpiVal  = 200
    legs = []
    for a in fig.get_axes():
        addLeg = a.get_legend()
        if not addLeg is None: legs.append(a.get_legend())
    if mode == 'png' or mode == 'all':
        fig.savefig(fileName+'.png',format='png',pad_inches=padding,bbox_inches='tight',
                    dpi=dpiVal,bbox_extra_artists=legs)
    if mode == 'eps': # or mode == 'all':
        fig.savefig(fileName+'.eps',format='eps',pad_inches=padding,bbox_inches='tight',
                    dpi=dpiVal,bbox_extra_artists=legs)
    if mode == 'pdf' or mode == 'all':
        fig.savefig(fileName+'.pdf',format='pdf',pad_inches=padding,bbox_inches='tight',
                    dpi=dpiVal,bbox_extra_artists=legs)
# ---------------------

def fig_format(ax, mode='save'):
    assert type(mode) == type('')
    assert mode.lower() in ['save','show'], 'Unknown mode'
    titleSize     = 24 #40 #38
    axLabelSize   = 20 #38 #36
    tickLabelSize = 18 #30 #28
    legendSize    = tickLabelSize-2
    textSize      = legendSize-2
    deltaShow     = 4

    def myformat(myax):
        if mode.lower() == 'show':
            for i in myax.get_children(): # Gets EVERYTHING!
                if isinstance(i, txt.Text):
                    i.set_size(textSize+3*deltaShow)

            for i in myax.get_lines():
                if i.get_marker() == 'D': continue # Don't modify baseline diamond
                i.set_linewidth(4)
                #i.set_markeredgewidth(4)
                i.set_markersize(10)

            leg = myax.get_legend()
            if not leg is None:
                for t in leg.get_texts(): t.set_fontsize(legendSize+deltaShow+6)
                th = leg.get_title()
                if not th is None:
                    th.set_fontsize(legendSize+deltaShow+6)

            myax.set_title(myax.get_title(),size=titleSize+deltaShow,weight='bold')
            myax.set_xlabel(myax.get_xlabel(),size=axLabelSize+deltaShow,weight='bold')
            myax.set_ylabel(myax.get_ylabel(),size=axLabelSize+deltaShow,weight='bold')
            myax.tick_params(labelsize=tickLabelSize+deltaShow)
            myax.patch.set_linewidth(3)
            for i in myax.get_xticklabels():
                i.set_size(tickLabelSize+deltaShow)
            for i in myax.get_xticklines():
                i.set_linewidth(3)
            for i in myax.get_yticklabels():
                i.set_size(tickLabelSize+deltaShow)
            for i in myax.get_yticklines():
                i.set_linewidth(3)

        elif mode.lower() == 'save':
            for i in myax.get_children(): # Gets EVERYTHING!
                if isinstance(i, txt.Text):
                    i.set_size(textSize)

            for i in myax.get_lines():
                if i.get_marker() == 'D': continue # Don't modify baseline diamond
                i.set_linewidth(4)
                #i.set_markeredgewidth(4)
                i.set_markersize(10)

            leg = myax.get_legend()
            if not leg is None:
                for t in leg.get_texts(): t.set_fontsize(legendSize)
                th = leg.get_title()
                if not th is None:
                    th.set_fontsize(legendSize)

            myax.set_title(myax.get_title(),size=titleSize,weight='bold')
            myax.set_xlabel(myax.get_xlabel(),size=axLabelSize,weight='bold')
            myax.set_ylabel(myax.get_ylabel(),size=axLabelSize,weight='bold')
            myax.tick_params(labelsize=tickLabelSize)
            myax.patch.set_linewidth(3)
            for i in myax.get_xticklabels():
                i.set_size(tickLabelSize)
            for i in myax.get_xticklines():
                i.set_linewidth(3)
            for i in myax.get_yticklabels():
                i.set_size(tickLabelSize)
            for i in myax.get_yticklines():
                i.set_linewidth(3)

    if type(ax) == type([]):
        for i in ax: myformat(i)
    else:
        myformat(ax)
        
# ---------------------
