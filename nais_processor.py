import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from datetime import datetime, timedelta
import matplotlib.dates as dts
import pandas as pd
import os
import locale
import warnings
import yaml
import re
from dateutil.parser import parse

dp_ion_lowres = np.array([7.949610066873187275e-01,
9.181737924552214603e-01,
1.060513600503926179e+00,
1.224959679823698799e+00,
1.414958699738506631e+00,
1.634499249798819331e+00,
1.888198514085806856e+00,
2.181403433339687226e+00,
2.520308747865528165e+00,
2.912095102815642989e+00,
3.365090891236600878e+00,
3.888962384293289887e+00,
4.494937535166431353e+00,
5.196070414640996837e+00,
6.007554438162747701e+00,
6.947095098447752193e+00,
8.035355151375323857e+00,
9.296489193192451594e+00,
1.075878902024538242e+01,
1.245546773082500103e+01,
1.442561898219513949e+01,
1.671539984850161886e+01,
1.937950186998520152e+01,
2.248299804137784363e+01,
2.610368545677439300e+01,
3.033508982931992648e+01,
3.529036394466827886e+01,
4.110740875515996606e+01])

dp_ion_hires = np.array([7.949610066875409942e-01,
8.241182580927455259e-01,
8.543462802115115995e-01,
8.856844966190522417e-01,
9.181737924559271180e-01,
9.518565575670346890e-01,
9.867767691571656119e-01,
1.022980023143976958e+00,
1.060513600503939280e+00,
1.099426538829467725e+00,
1.139769706974383290e+00,
1.181595851620098836e+00,
1.224959679823832692e+00,
1.269917926237073669e+00,
1.316529439401527446e+00,
1.364855258822949002e+00,
1.414958699738505521e+00,
1.466905426169209825e+00,
1.520763551055280605e+00,
1.576603729908572893e+00,
1.634499249758283312e+00,
1.694526135045845150e+00,
1.756763254297374122e+00,
1.821292427408487402e+00,
1.888198514077982448e+00,
1.957569584416486874e+00,
2.029496978295826093e+00,
2.104075485232879572e+00,
2.181403433339692999e+00,
2.261582859332240680e+00,
2.344719627899444880e+00,
2.430923596381127538e+00,
2.520308748234819429e+00,
2.612993381959131334e+00,
2.709100248607977601e+00,
2.808756747008758214e+00,
2.912095102815391190e+00,
3.019252527502206185e+00,
3.130371476488065241e+00,
3.245599771518457466e+00,
3.365090891236902415e+00,
3.489004141625421163e+00,
3.617504913802647604e+00,
3.750764898980823325e+00,
3.888962384292350194e+00,
4.032282462010540414e+00,
4.180917334988799361e+00,
4.335066596457491706e+00,
4.494937535166420695e+00,
4.660745425685828280e+00,
4.832713857244105959e+00,
5.011075095921604827e+00,
5.196070414634600176e+00,
5.387950447963188338e+00,
5.586975636988666061e+00,
5.793416574994243007e+00,
6.007554436026940614e+00,
6.229681456686657626e+00,
6.460101359513265251e+00,
6.699129841480823799e+00,
6.947095098447658934e+00,
7.204338334753603412e+00,
7.471214336832219693e+00,
7.748092044569639292e+00,
8.035355151259254924e+00,
8.333402769720130721e+00,
8.642650082548662738e+00,
8.963529060188569986e+00,
9.296489193192432055e+00,
9.641998264711249433e+00,
1.000054318925402619e+01,
1.037263084818457415e+01,
1.075878902020633987e+01,
1.115956729412175541e+01,
1.157553811904305796e+01,
1.200729781901941351e+01,
1.245546773067095359e+01,
1.292069536930353912e+01,
1.340365566604162062e+01,
1.390505228292856543e+01,
1.442561898541390164e+01,
1.496612109930933521e+01,
1.552735707366413109e+01,
1.611016009905102919e+01,
1.671539984589247041e+01,
1.734398428907315903e+01,
1.799686165161790541e+01,
1.867502249416057936e+01,
1.937950186998513757e+01,
2.011138166759334922e+01,
2.087179308880621065e+01,
2.166191924893708176e+01,
2.248299804205197461e+01,
2.333632504917748918e+01,
2.422325676528101823e+01,
2.514521399524151590e+01,
2.610368545923754624e+01,
2.710023166043850651e+01,
2.813648906294003282e+01,
2.921417447477816509e+01,
3.033508982951557797e+01,
3.150112726419725817e+01,
3.271427458287274703e+01,
3.397662110330684015e+01,
3.529036394400620225e+01,
3.665781484216577724e+01,
3.808140737765210559e+01,
3.956370485580886509e+01,
4.110740875118314364e+01])

dp_par_lowres = np.array([7.498942093324539870e-01,
8.659643233600640144e-01,
9.999999999999980016e-01,
1.154781984689456031e+00,
1.333521432163321974e+00,
1.539926526059490097e+00,
1.778279410038920094e+00,
2.053525026457140079e+00,
2.371373705661659947e+00,
2.738419634264360081e+00,
3.162277660168379967e+00,
3.651741272548380213e+00,
4.216965034285819591e+00,
4.869675251658620141e+00,
5.623413251903479626e+00,
6.493816315762099833e+00,
7.498942093324560076e+00,
8.659643233600640144e+00,
1.000000000000000000e+01,
1.154781984689457985e+01,
1.333521432163323972e+01,
1.539926526059490008e+01,
1.778279410038922137e+01,
2.053525026457139901e+01,
2.371373705661660125e+01,
2.738419634264360170e+01,
3.162277660168379967e+01,
3.651741272548380124e+01,
4.216965034285819769e+01])

dp_par_hires = np.array([7.498942093324509894e-01,
7.773650302387707933e-01,
8.058421877614766471e-01,
8.353625469578208618e-01,
8.659643233600599066e-01,
8.976871324473085778e-01,
9.305720409296931450e-01,
9.646616199111932577e-01,
9.999999999999937828e-01,
1.036632928437691614e+00,
1.074607828321310965e+00,
1.113973859994795701e+00,
1.154781984689451368e+00,
1.197085030495722791e+00,
1.240937760751712249e+00,
1.286396944936966991e+00,
1.333521432163316423e+00,
1.382372227357891781e+00,
1.433012570236954719e+00,
1.485508017172767037e+00,
1.539926526059483658e+00,
1.596338544287933647e+00,
1.654817099943172609e+00,
1.715437896342869806e+00,
1.778279410038913433e+00,
1.843422992409100791e+00,
1.910952974970430596e+00,
1.980956778550328590e+00,
2.053525026457135638e+00,
2.128751661796361994e+00,
2.206734069084578920e+00,
2.287573200318384625e+00,
2.371373705661643960e+00,
2.458244068920186098e+00,
2.548296747979334587e+00,
2.641648320386080329e+00,
2.738419634264348979e+00,
2.838735964758742458e+00,
2.942727176209269491e+00,
3.050527890267013209e+00,
3.162277660168366644e+00,
3.278121151393445398e+00,
3.398208328942545542e+00,
3.522694651473087468e+00,
3.651741272548362449e+00,
3.785515249258615267e+00,
3.924189758484520674e+00,
4.067944321083031944e+00,
4.216965034285806269e+00,
4.371444812611072983e+00,
4.531583637600800962e+00,
4.697588816706474546e+00,
4.869675251658613035e+00,
5.048065716667452740e+00,
5.232991146814928385e+00,
5.424690937011307668e+00,
5.623413251903471632e+00,
5.829415347136054137e+00,
6.042963902381307761e+00,
6.264335366568834829e+00,
6.493816315762091840e+00,
6.731703824144960713e+00,
6.978305848598641781e+00,
7.233941627366726301e+00,
7.498942093324536096e+00,
7.773650302387736133e+00,
8.058421877614796003e+00,
8.353625469578238594e+00,
8.659643233600629486e+00,
8.976871324473117753e+00,
9.305720409296965201e+00,
9.646616199111967660e+00,
9.999999999999975131e+00,
1.036632928437695433e+01,
1.074607828321314962e+01,
1.113973859994799831e+01,
1.154781984689455676e+01,
1.197085030495727231e+01,
1.240937760751716823e+01,
1.286396944936971742e+01,
1.333521432163321130e+01,
1.382372227357896577e+01,
1.433012570236959782e+01,
1.485508017172772099e+01,
1.539926526059488943e+01,
1.596338544287939420e+01,
1.654817099943178249e+01,
1.715437896342875845e+01,
1.778279410038919650e+01,
1.843422992409107408e+01,
1.910952974970437523e+01,
1.980956778550335784e+01,
2.053525026457143454e+01,
2.128751661796369632e+01,
2.206734069084587091e+01,
2.287573200318393063e+01,
2.371373705661652664e+01,
2.458244068920194891e+01,
2.548296747979344090e+01,
2.641648320386089921e+01,
2.738419634264358749e+01,
2.838735964758752317e+01,
2.942727176209279705e+01,
3.050527890267023423e+01,
3.162277660168377835e+01,
3.278121151393457211e+01,
3.398208328942558154e+01,
3.522694651473100436e+01,
3.651741272548376571e+01,
3.785515249258629922e+01,
3.924189758484535417e+01,
4.067944321083047043e+01,
4.216965034285821901e+01])

mob_ion_lowres = np.array([3.162277660168379937e-04,
2.371373705661659990e-04,
1.778279410038920258e-04,
1.333521432163320159e-04,
1.000000000000000048e-04,
7.498942093324559917e-05,
5.623413251903490022e-05,
4.216965034285820205e-05,
3.162277660168380208e-05,
2.371373705661660125e-05,
1.778279410038919852e-05,
1.333521432163319990e-05,
1.000000000000000082e-05,
7.498942093324561442e-06,
5.623413251903490361e-06,
4.216965034285830030e-06,
3.162277660168380038e-06,
2.371373705661659871e-06,
1.778279410038920148e-06,
1.333521432163330027e-06,
1.000000000000000167e-06,
7.498942093324570124e-07,
5.623413251903499890e-07,
4.216965034285829924e-07,
3.162277660168379721e-07,
2.371373705661660136e-07,
1.778279410038920042e-07,
1.333521432163329868e-07])

mob_ion_hires = np.array([3.162277660168379394e-04,
2.942727176209282208e-04,
2.738419634264362135e-04,
2.548296747979347430e-04,
2.371373705661656195e-04,
2.206734069084591093e-04,
2.053525026457147184e-04,
1.910952974970441833e-04,
1.778279410038924053e-04,
1.654817099943182849e-04,
1.539926526059493548e-04,
1.433012570236964366e-04,
1.333521432163325580e-04,
1.240937760751721229e-04,
1.154781984689459841e-04,
1.074607828321319270e-04,
1.000000000000001539e-04,
9.305720409297006353e-05,
8.659643233600669858e-05,
8.058421877614834719e-05,
7.498942093324574825e-05,
6.978305848598679792e-05,
6.493816315762129137e-05,
6.042963902381343169e-05,
5.623413251903505608e-05,
5.232991146814961977e-05,
4.869675251658644782e-05,
4.531583637600830909e-05,
4.216965034285835112e-05,
3.924189758484547866e-05,
3.651741272548388522e-05,
3.398208328942570449e-05,
3.162277660168389694e-05,
2.942727176209291627e-05,
2.738419634264370944e-05,
2.548296747979355561e-05,
2.371373705661663852e-05,
2.206734069084598140e-05,
2.053525026457154028e-05,
1.910952974970448067e-05,
1.778279410038930016e-05,
1.654817099943188405e-05,
1.539926526059498630e-05,
1.433012570236968973e-05,
1.333521432163329985e-05,
1.240937760751725159e-05,
1.154781984689463535e-05,
1.074607828321322658e-05,
1.000000000000004825e-05,
9.305720409297036169e-06,
8.659643233600697641e-06,
8.058421877614859791e-06,
7.498942093324598711e-06,
6.978305848598699613e-06,
6.493816315762147941e-06,
6.042963902381361465e-06,
5.623413251903522549e-06,
5.232991146814976546e-06,
4.869675251658659012e-06,
4.531583637600844292e-06,
4.216965034285847818e-06,
3.924189758484559555e-06,
3.651741272548399110e-06,
3.398208328942580868e-06,
3.162277660168399520e-06,
2.942727176209300267e-06,
2.738419634264378652e-06,
2.548296747979362761e-06,
2.371373705661670459e-06,
2.206734069084603985e-06,
2.053525026457159618e-06,
1.910952974970453149e-06,
1.778279410038934759e-06,
1.654817099943192683e-06,
1.539926526059502526e-06,
1.433012570236972658e-06,
1.333521432163333415e-06,
1.240937760751728420e-06,
1.154781984689466414e-06,
1.074607828321325199e-06,
1.000000000000007155e-06,
9.305720409297056286e-07,
8.659643233600717122e-07,
8.058421877614878214e-07,
7.498942093324614593e-07,
6.978305848598716341e-07,
6.493816315762163400e-07,
6.042963902381375018e-07,
5.623413251903534830e-07,
5.232991146814988828e-07,
4.869675251658670447e-07,
4.531583637600854986e-07,
4.216965034285857982e-07,
3.924189758484568661e-07,
3.651741272548407898e-07,
3.398208328942588279e-07,
3.162277660168406720e-07,
2.942727176209307466e-07,
2.738419634264385216e-07,
2.548296747979369008e-07,
2.371373705661676282e-07,
2.206734069084609914e-07,
2.053525026457164753e-07,
1.910952974970458125e-07,
1.778279410038939365e-07,
1.654817099943196971e-07,
1.539926526059506708e-07,
1.433012570236976523e-07,
1.333521432163337015e-07])

def visc(temp):
    """ Calculate viscosity of air """

    nyy_ref=18.203e-6
    S=110.4
    temp_ref=293.15
    nyy=nyy_ref*((temp_ref+S)/(temp+S))*((temp/temp_ref)**(3./2.))
    return nyy

def rlambda(temp,press):
    """ Calculate mean-free path """

    kaasuv=8.3143
    dm=3.7e-10
    avoc=6.022e23
    return kaasuv*temp/(np.sqrt(2.)*avoc*press*np.pi*dm*dm)

def cunn(Dp,temp,press):
    """ Calculate Cunningham correction """

    lambd = rlambda(temp,press)
    return 1. + 2.*lambd/Dp * (1.165 + 0.483 * np.exp(-0.997*Dp/lambd))

def diffuus(dpp,temp,press):
    """ Calculate diffusion coefficient """

    K=1.38e-23
    return (K*temp*cunn(dpp,temp,press))/(3.*np.pi*visc(temp)*dpp)

def tubeloss(dpp,pflow,plength,temp,press):
    """ Laminar diffusion losses in circular straight conduit """

    rmuu = np.pi*diffuus(dpp,temp,press)*plength/pflow;
    pene = np.zeros(rmuu.shape)
    for i in range(0,len(dpp)):
        if rmuu[i] < 0.02:
            pene[i]=1. - 2.56*rmuu[i]**(2./3.) \
                    + 1.2*rmuu[i]+0.177*rmuu[i]**(4./3.)
        else:
            pene[i]=0.819*np.exp(-3.657*rmuu[i])+0.097 \
                    *np.exp(-22.3*rmuu[i])+0.032*np.exp(-57.*rmuu[i])
    return pene

def x2dlogx(x):
    """ Calculate log-differences for monotonically
    increasing or decreasing vector x """

    logx = np.log10(x)
    logx_mid = (logx[1:] + logx[:-1])/2.0 
    logx_mid_first_value = logx[0] + (logx[0] - logx_mid[0])
    logx_mid_last_value  = logx[-1] - (logx_mid[-1] - logx[-1])
    logx_mid = np.insert(logx_mid,0,logx_mid_first_value)
    logx_mid = np.append(logx_mid,logx_mid_last_value)
    dlogx = np.abs(np.diff(logx_mid))
    return dlogx

def datetime2datenum(dt):
    """ Convert from python datetime to matlab datenum """

    ord = dt.toordinal()
    mdn = dt + timedelta(days = 366)
    frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0,tzinfo=dt.tzinfo)).seconds \
           / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac

def datenum2datetime(matlab_datenum):
    """ Convert from matlab datenum to python datetime """

    return datetime.fromordinal(int(matlab_datenum)) \
    + timedelta(days=matlab_datenum%1) \
    - timedelta(days = 366)


def plot_sumfile(handle,v,clim=(10,100000)):
    """ Plot UHEL's sum-formatted aerosol number-size distribution """
    
    plt.rcParams.update({'font.size': 10,
                         'axes.titlesize': 10,
                         'axes.labelsize': 10,
                         'xtick.labelsize': 10,
                         'ytick.labelsize': 10,
                         'figure.titlesize': 10,
                         'legend.fontsize': 10})
    
    time = v[1:,0]
    dp = v[0,2:]
    data = v[1:,2:]
    mesh_dp,mesh_time = np.meshgrid(dp,time)
    pcolorplot = handle.pcolormesh(mesh_time,mesh_dp,data,
                                   norm=colors.LogNorm(),
    linewidth=0,rasterized=True)
    handle.set_yscale('log')
    pcolorplot.set_clim(clim)
    pcolorplot.set_edgecolor('face')
    handle.autoscale(tight='true')

    handle.grid('on',which='both',linestyle='--',color='w',lw=0.5)
    handle.xaxis.set_major_locator(dts.HourLocator(interval=1))
    handle.xaxis.set_major_formatter(dts.DateFormatter('%H'))
    plt.setp(handle.get_xticklabels(), rotation=80)
    handle.set_ylabel('Dp, [m]')
    handle.set_xlabel('Computer time, [hours]')
    cbar = plt.colorbar(pcolorplot,ax=handle,
                        ticks=LogLocator(subs=range(10)))
    cbar.set_label('dN/dlogDp, [cm-3]')
    return pcolorplot

def nais_processor(load_path = '/path/to/raw/data/',
                   save_path = '/path/to/processed/files/',
                   config_file = '/path/to/config_file.yml',
                   log_file = '/path/to/log_file.log'):
    """ Function that processes data from the NAIS 

    The function reads the raw data files from Neutral cluster and 
    Air Ion Spectrometer (NAIS), applies corrections (diffusion losses in 
    the inlet line, conversion to standard conditions and R. Wagner's 
    ion mode calibration) to the measured number concentrations and 
    saves the data as a University of Helsinki sum-formatted number-size
    distribution.

    When the function is called it will process all the unprocessed 
    raw files in the source folder.

    Sum-format:
        [1:,0] = Time (MATLAB datenum) 
        [0,2:] = Geometric mean diameter of size-channel (m)
        [1:,1] = Integrated total number concentration (cm-3)
        [1:,2:] = Normalized number concentrations (cm-3)

    Example configuration file:
        # NAIS-5-1.yml

        model: NAIS-5-1
        inverter_resolution: 'high'
        sampleflow: 54.0 # L/min
        date_format: '%Y%m%d'
        particle_files: '-block-particles.spectra'
        ion_files: '-block-ions.spectra'
        diagnostic_files: '-block.records'
        delimiter: '\t'
        temperature_name: 'temperature.mean'
        pressure_name: 'baro.mean'
        pipe_length: 0.5  # meters

    Arguments:
        load_path (str): path from where data is loaded
        save_path (str): path where data is saved
        config_file (str): path to configuration file
    
    Example use:
        import sys
        sys.path.append('/home/user/nais_processor')
        from nais_processor import *
        nais_processor(load_path='/home/user/NAIS/data/',
                           save_path='/home/user/NAIS/processed/',
                           config_file='/home/user/NAIS/NAIS-5-1.yml',
                           log_file='/home/user/NAIS/NAIS-5-1.log')
        nais_plotter(load_path='/home/user/NAIS/processed/',
                           save_path='/home/user/NAIS/figs/')
    """
    
    # Ignore all warnings
    warnings.filterwarnings("ignore")

    # Try to read the configuration file
    with open(config_file,'r') as stream:
        try:
            print 'config file: '+config_file
            config = yaml.load(stream)
            pressure_name = config['pressure_name']
            temperature_name = config['temperature_name']
            inverter_reso = config['inverter_resolution']
            sampleflow = config['sampleflow']
            date_format = config['date_format']
            particle_files = config['particle_files']
            ion_files = config['ion_files']
            diagnostic_files = config['diagnostic_files']
            delimiter = config['delimiter']
            model = config['model']
            pipelength = config['pipe_length']
        except Exception as error:
            print 'error: '+str(error)
            return

    # List files in the load and save directories
    load_filenames = os.listdir(load_path)
    save_filenames = os.listdir(save_path)

    # Find the dates that need to be processed
    load_dates = np.array([datetime.strptime(re.search('.+?(?='+ion_files+')',x).group(0),\
                 date_format) for x in load_filenames if ion_files in x])
    save_dates = np.array([datetime.strptime(re.search('(?<='+model+'n)(.*)(?=nds.sum)',x).group(0),\
                 '%Y%m%d') for x in save_filenames if re.match(model+'n........nds.sum',x)])
    datetimes = np.sort(np.setdiff1d(load_dates,save_dates))
    date_strings = np.array([datetime.strftime(x,date_format) for x in datetimes])

    # Check if any new data files found
    if len(date_strings)==0:
        return
    else:
        # Do not process the latest file if its from today
        today = datetime.today().strftime(date_format)
        last_day = date_strings[-1]
        if last_day==today:
            last_index = len(datetimes)-1
        else:
            last_index = len(datetimes)

    # Define standard conditions
    temp_ref = 273.15 # K
    pres_ref = 101325.0 # Pa

    # Determine the inverter resolution
    if inverter_reso=='low':
         dp_ion = dp_ion_lowres
         mob_ion = mob_ion_lowres
         dp_par = dp_par_lowres
    elif inverter_reso=='high':
         dp_ion = dp_ion_hires
         mob_ion = mob_ion_hires
         dp_par = dp_par_hires
    else:
        print "error: 'inverter resolution must be 'low' or 'high''"
        return

    # Calculate the log-differences
    dlogdp_ion = x2dlogx(dp_ion)
    dlogmob_ion = x2dlogx(mob_ion)
    dlogdp_par = x2dlogx(dp_par)

    # Run through the unprocessed dates
    for i in range(0,last_index):

        # If model is NAIS then particles are measured
        if model[:4]=='NAIS':
      
            try:
                print 'processing '+date_strings[i]
            
                # Read the ion data
                ions = pd.read_table(load_path + date_strings[i] + ion_files,
                                     sep=delimiter,
                                     comment='#',
                                     engine='python',
                                     error_bad_lines=False,
                                     header=None)

                # Read the particle data 
                particles = pd.read_table(load_path + date_strings[i] + particle_files,
                                          sep=delimiter,
                                          comment='#',
                                          engine='python',
                                          error_bad_lines=False,
                                          header=None)

                # Read the diagnostic data
                records = pd.read_table(load_path + date_strings[i] + diagnostic_files,
                                        sep=delimiter,
                                        comment='#',
                                        engine='python',
                                        error_bad_lines=False,
                                        header=None)

                # Remove duplicate rows based on first column (time)
                ions = ions.drop_duplicates(subset=None,
                                            keep='first',
                                            inplace=False)

                particles = particles.drop_duplicates(subset = None,
                                                      keep='first',
                                                      inplace=False)

                records = records.drop_duplicates(subset=None,
                                                  keep='first',
                                                  inplace=False)

                # Set the first row as the header and remove it from the actual data
                ions.columns = ions.iloc[0,:]
                ions = ions.reindex(ions.index.drop(0))

                particles.columns = particles.iloc[0,:]
                particles = particles.reindex(particles.index.drop(0))

                records.columns = records.iloc[0,:]
                records = records.reindex(records.index.drop(0))

                # Calculate the ion number-size distribution
                ions = ions.set_index(ions.columns[0])
                ions.index = [parse(x) for x in ions.index]
                neg_ions = ions.iloc[:,3:3+len(dp_ion)].astype(float).interpolate().values
                pos_ions = ions.iloc[:,3+2*len(dp_ion):3+3*len(dp_ion)].astype(float).interpolate().values
                neg_ions = neg_ions * dlogmob_ion / dlogdp_ion
                pos_ions = pos_ions * dlogmob_ion / dlogdp_ion

                # Calculate the particle number-size distribution
                particles = particles.set_index(particles.columns[0])
                particles.index = [parse(x) for x in particles.index]
                neg_particles = particles.iloc[:,3:3+len(dp_par)].astype(float).interpolate().values
                pos_particles = particles.iloc[:,3+2*len(dp_par):3+3*len(dp_par)].astype(float).interpolate().values

                # Index records by the operation mode
                records = records.set_index('opmode')

                # Then extract the ion records
                ion_records = records.loc['ions'].set_index(records.columns[0])
                ion_records.index = [parse(x) for x in ion_records.index]
                ion_records = ion_records.reindex(index=ions.index)

                # And the particle records
                particle_records = records.loc['particles'].set_index(records.columns[0])
                particle_records.index = [parse(x) for x in particle_records.index]
                particle_records = particle_records.reindex(index=particles.index)

                # List all variable names in the diagnostic file
                diag_params = list(records)

                # Extract temperature data if it exists, else use standard conditions
                if temperature_name=='':
                    temp_ions = temp_ref*np.ones(neg_ions.shape[0])
                    temp_particles = temp_ref*np.ones(neg_particles.shape[0])
                elif temperature_name in diag_params:
                    temp_ions = 273.15 + ion_records[temperature_name].astype(float).interpolate().values.flatten()
                    temp_particles = 273.15 + particle_records[temperature_name].astype(float).interpolate().values.flatten()
                else:
                    err_msg = 'name '+"'"+temperature_name+"'"+' not found in '+"'"+date_strings[i]+diagnostic_files+"'"
                    print 'error: '+err_msg
                    with open(log_file,'a+') as the_log_file:
                        the_log_file.write('%s\nerror while processing: %s\n%s\n' % (datetime.today().strftime('%Y-%m-%d %H:%M'),date_strings[i],err_msg))
                    continue
          
                # Extract pressure data if it exists, else use standard conditions
                if pressure_name=='':
                    pres_ions = pres_ref*np.ones(neg_ions.shape[0])
                    pres_particles = pres_ref*np.ones(neg_particles.shape[0])
                elif pressure_name in diag_params:
                    pres_ions = 100.0*ion_records[pressure_name].astype(float).interpolate().values.flatten()
                    pres_particles = 100.0*particle_records[pressure_name].astype(float).interpolate().values.flatten()
                else:
                    err_msg = 'name '+"'"+pressure_name+"'"+' not found in '+"'"+date_strings[i]+diagnostic_files+"'"
                    print 'error: '+err_msg
                    with open(log_file,'a+') as the_log_file:
                        the_log_file.write('%s\nerror while processing: %s\n%s\n' % (datetime.today().strftime('%Y-%m-%d %H:%M'),date_strings[i],err_msg))
                    continue

                # Pressure and temperature correction
                stp_corr_ions = (pres_ref*temp_ions)/(temp_ref*pres_ions)
                neg_ions = (stp_corr_ions*neg_ions.T).T
                pos_ions = (stp_corr_ions*pos_ions.T).T

                stp_corr_particles = (pres_ref*temp_particles)/(temp_ref*pres_particles)
                neg_particles = (stp_corr_particles*neg_particles.T).T
                pos_particles = (stp_corr_particles*pos_particles.T).T

                # Diffusion loss correction
                throughput_ions = np.zeros(neg_ions.shape)
                for j in range(0,throughput_ions.shape[0]):
                    throughput_ions[j,:] = tubeloss(dp_ion*1e-9,
                                                    sampleflow*1.667e-5,
                                                    pipelength,
                                                    temp_ions[j],
                                                    pres_ions[j])

                throughput_particles = np.zeros(neg_particles.shape)
                for j in range(0,throughput_particles.shape[0]):
                    throughput_particles[j,:] = tubeloss(dp_par*1e-9,
                                                         sampleflow*1.667e-5,
                                                         pipelength,
                                                         temp_particles[j],
                                                         pres_particles[j])

                neg_ions = throughput_ions*neg_ions
                pos_ions = throughput_ions*pos_ions
                neg_particles = throughput_particles*neg_particles
                pos_particles = throughput_particles*pos_particles

                # Robert Wagner's calibration (only ions)
                roberts_corr = 0.713*dp_ion**0.120
                neg_ions = roberts_corr*neg_ions
                pos_ions = roberts_corr*pos_ions

                # Integrate total number concentrations
                total_neg_ions = np.nansum(neg_ions*dlogdp_ion,axis=1)[np.newaxis].T
                total_pos_ions = np.nansum(pos_ions*dlogdp_ion,axis=1)[np.newaxis].T
                total_neg_particles = np.nansum(neg_particles*dlogdp_par,axis=1)[np.newaxis].T
                total_pos_particles = np.nansum(pos_particles*dlogdp_par,axis=1)[np.newaxis].T

                # Define time vectors
                time_ions = np.array([datetime2datenum(x) for x in ions.index])[np.newaxis].T
                time_particles = np.array([datetime2datenum(x) for x in particles.index])[np.newaxis].T
                
                # Construct the headers
                ion_header = np.insert(dp_ion*1e-9,(0,0),0)[np.newaxis]
                particle_header = np.insert(dp_par*1e-9,(0,0),0)[np.newaxis]

                # Construct the sum-files
                negions = np.concatenate((ion_header,np.concatenate((time_ions,total_neg_ions,neg_ions),axis=1)))
                posions = np.concatenate((ion_header,np.concatenate((time_ions,total_pos_ions,pos_ions),axis=1)))
                negparticles = np.concatenate((particle_header,np.concatenate\
                                                    ((time_particles,total_neg_particles,neg_particles),axis=1)))
                posparticles = np.concatenate((particle_header,np.concatenate\
                                                    ((time_particles,total_pos_particles,pos_particles),axis=1)))

                # Save the sum matrices using the standard names
                ds = datetime.strftime(datetimes[i],'%Y%m%d')
                np.savetxt(save_path+model+'n'+ds+'nds.sum',negions)
                np.savetxt(save_path+model+'p'+ds+'nds.sum',posions)
                np.savetxt(save_path+model+'n'+ds+'np.sum',negparticles)
                np.savetxt(save_path+model+'p'+ds+'np.sum',posparticles)

            except Exception as error:
                print 'error: '+str(error)
                with open(log_file,'a+') as the_log_file:
                    the_log_file.write("%s\nunexpected error while processing: %s\n%s\n" % (datetime.today().strftime('%Y-%m-%d %H:%M'),date_strings[i],error))
                continue

        else: 
            continue

def nais_plotter(load_path='/path/to/processed/files/',
                 save_path='/path/to/figures/',
                 config_file='/path/to/config_file.yml'):
    """ Function to plot the processed NAIS data """

    # Try to read the configuration file
    with open(config_file,'r') as stream:
        try:
            config = yaml.load(stream)
            model = config['model']
            location = config['location']
        except Exception as error:
            print 'error: '+str(error)
            return

    plt.style.use('dark_background')

    # List files in the load and save directories
    load_filenames = os.listdir(load_path)
    save_filenames = os.listdir(save_path)
    
    # Find the dates that need to be plotted
    load_dates = np.array([datetime.strptime(re.search('(?<='+model+'n)(.*)(?=nds.sum)',x).group(0),\
                 '%Y%m%d') for x in load_filenames if re.match(model+'n........nds.sum',x)])
    save_dates = np.array([datetime.strptime(re.search('(?<='+model+'_)(.*)(?=.png)',x).group(0),\
                 '%Y%m%d') for x in save_filenames if re.match(model+'_.........png',x)])
    datetimes = np.setdiff1d(load_dates,save_dates)
    date_strings = [datetime.strftime(x,'%Y%m%d') for x in datetimes]

    for i in range(0,len(datetimes)):

        # If model is NAIS then particles are measured
        if model[:4]=='NAIS':
            
            print "plotting " + date_strings[i]
            
            pos_ion_data = np.loadtxt(load_path+model+'p'+date_strings[i]+'nds.sum')
            neg_ion_data = np.loadtxt(load_path+model+'n'+date_strings[i]+'nds.sum')
            pos_particle_data = np.loadtxt(load_path+model+'p'+date_strings[i]+'np.sum')
            neg_particle_data = np.loadtxt(load_path+model+'n'+date_strings[i]+'np.sum')

            fig,ax = plt.subplots(2,2,figsize=(9,6))
            ax = ax.flatten()

            plot_sumfile(ax[0],pos_ion_data,clim=(10,10000))
            plot_sumfile(ax[1],neg_ion_data,clim=(10,10000))
            plot_sumfile(ax[2],pos_particle_data,clim=(10,100000))
            plot_sumfile(ax[3],neg_particle_data,clim=(10,100000))
                       
            ax[0].set_title('Positive ions')
            ax[1].set_title('Negative ions')
            ax[2].set_title('Particles (positive polarity)')
            ax[3].set_title('Particles (negative polarity)')

            plt.tight_layout()

            fig.suptitle(model + ' ' + datetimes[i].strftime('%Y-%m-%d') + ' ' + location, y=1)

            plt.savefig(save_path+model+'_'+date_strings[i]+'.png',
                        bbox_inches='tight',dpi=300)

            plt.close()

        else:
            continue
