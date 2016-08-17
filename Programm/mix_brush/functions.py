# -*- coding: utf-8 -*-

import numpy as np

#############################
### Some useful functions ###
#############################

#~ A derivative by 7 datapoints
def dervec7(v,h=1):
    global f
    # Находит производную по 7и точкам
    # v - функция, призводную которой ищем
    # h - шаг
    n=len(v);
    f=zeros(n)
    for k in range(n):
        if k==0:
            f[0]=(-147.0*v[0]+360.0*v[1]-450.0*v[2]+400.0*v[3]-225.0*v[4]+72.0*v[5]-10.0*v[6])/60.0/h ;
        elif k==1:
            f[1]=(-10.0*v[0]-77.0*v[1]+150.0*v[2]-100.0*v[3]+50.0*v[4]-15.0*v[5]+2.0*v[6])/60.0/h ;
        elif k==2:
            f[2]=(2.0*v[0]-24.0*v[1]-35.0*v[2]+80.0*v[3]-30.0*v[4]+8.0*v[5]-v[6])/60.0/h ;
        elif k > 2 and k < n-3:
            f[k]=(-v[k-3]+9.0*v[k-2]-45.0*v[k-1]+45.0*v[k+1]-9.0*v[k+2]+v[k+3])/60.0/h;
        elif k==n-3:
            f[k]=(v[k-4]-8.0*v[k-3]+30.0*v[k-2]-80.0*v[k-1]+35.0*v[k]+24.0*v[k+1]-2.0*v[k+2])/60.0/h;
        elif k==n-2:
            f[k]=(-2.0*v[k-5]+15.0*v[k-4]-50.0*v[k-3]+100.0*v[k-2]-150.0*v[k-1]+77.0*v[k]+10.0*v[k+1])/60.0/h;
        else:
            f[k]=(10.0*v[k-6]-72.0*v[k-5]+225.0*v[k-4]-400.0*v[k-3]+450.0*v[k-2]-360.0*v[k-1]+147.0*v[k])/60.0/h;
    return f

#~ A derivative by 4 datapoints
def dervec4(v,h=1):
    global f,k
    # Находит производную по 4ем точкам
    # v - функция, призводную которой ищем
    # h - шаг
    n=len(v);
    f=zeros(n)
    for k in range(n):
        if k==0:
            f[0]=(-11.0*v[0]+18.0*v[1]-9.0*v[2]+2.0*v[3])/6.0/h ;
        elif k==1:
            f[1]=(-2.0*v[0]-3.0*v[1]+6.0*v[2]-v[3])/6.0/h ;
        elif k>1 and k<n-1:
            f[k]=(v[k-2]-6.0*v[k-1]+3.0*v[k]+2.0*v[k+1])/6.0/h ;
        else:
            f[k]=(-2*v[k-3]+9.0*v[k-2]-18.0*v[k-1]+11.0*v[k])/6.0/h;
    return f

#~ A derivative by 3 datapoints
def dervec3(v,h=1):
    global f
    # Находит производную по 3м точкам
    # v - функция, призводную которой ищем
    # h - шаг
    n=len(v);
    f=zeros(n)
    for k in range(n):
        if k==0:
            f[0]=(-3.0*v[0]+4.0*v[1]-v[2])/2.0/h ;
        elif k>0 and k<n-1:
            f[k]=(v[k+1]-v[k-1])/2.0/h ;
        else:
            f[k]=(v[k-2]-4.0*v[k-1]+3.0*v[k])/2.0/h ;
    return f

# ?? finds the first fixed point
def findmin(F, Rstar_range):
    # ruturns the position of minimum of F as function of R
    DF=dervec3(F)
    if sum(sign(DF)) == -len(DF):
        R0 = Rstar_range[-1]
    elif sum(sign(DF)) == len(DF):
        R0 = Rstar_range[0]
    else:
        # ?? do not understand what find() is doing
        nearest_to_zero_positive_index = find(sign(DF) ==  1)[ 0]
        nearest_to_zero_negative_index = find(sign(DF) == -1)[-1]
        R1 = Rstar_range[nearest_to_zero_positive_index];R1=float(R1)
        R2 = Rstar_range[nearest_to_zero_negative_index];R2=float(R2)
        DF1 = DF[nearest_to_zero_positive_index];DF1=float(DF1)
        DF2 = DF[nearest_to_zero_negative_index];DF2=float(DF2)
        R0 = -DF1/(DF2-DF1)*(R2-R1)+R1
    return R0
def findmin7(F, Rstar_range):
    # ruturns the position of minimum of F as function of R
    DF=dervec7(F)
    if sum(sign(DF)) == -len(DF):
        R0 = Rstar_range[-1]
    elif sum(sign(DF)) == len(DF):
        R0 = Rstar_range[0]
    else:
        i = 0
        while sign(DF[i]) ==  -1:
            index = i
            i = i+1
        nearest_to_zero_positive_index = i
        nearest_to_zero_negative_index = i-1
        R1 = Rstar_range[nearest_to_zero_positive_index];R1=float(R1)
        R2 = Rstar_range[nearest_to_zero_negative_index];R2=float(R2)
        DF1 = DF[nearest_to_zero_positive_index];DF1=float(DF1)
        DF2 = DF[nearest_to_zero_negative_index];DF2=float(DF2)
        R0 = -DF1/(DF2-DF1)*(R2-R1)+R1
    return R0



#~ Makes xy plot, xname and yname are the names of axis
def vplot(x,y,xname = 'x', yname = 'y', xlog = False, ylog = False, color = 'black', PlotLine = True, marker = 'circle', markersize = '2pt'):

    global g, graph, xy
    try:
        type (g) == veusz.Embedded
        try:
            type(graph) == veusz.WidgetNode
        except NameError:
            page = g.Root.page1
            graph = page.graph1
    except NameError:
        g = veusz.Embedded('F')
        g.EnableToolbar()
        page = g.Root.Add('page')
        graph = page.Add('graph')


    x_dataname  = xname
    y_dataname  = yname

    if len(np.shape(x)) == 2:
        x_data = x[0]
        x_data_err = x[1]
        g.SetData(x_dataname, x_data, symerr = x_data_err)
    else:
        x_data = x
        g.SetData(x_dataname, x_data)
    if len(np.shape(y)) == 2:
        y_data = y[0]
        y_data_err = y[1]
        g.SetData(y_dataname, y_data, symerr = y_data_err)
    else:
        y_data = y
        g.SetData(y_dataname, y_data)



    xy = graph.Add('xy')
    xy.xData.val = x_dataname
    xy.yData.val = y_dataname
    xy.marker.val = marker
    xy.MarkerFill.color.val = color
    xy.markerSize.val = markersize
    #~ xy_sfbox.ErrorBarLine.width.val = '2pt'
    #~ xy_sfbox.ErrorBarLine.transparency.val = 50
    xy.PlotLine.width.val = '2pt'
    xy.PlotLine.style.val = 'solid'
    xy.PlotLine.color.val = color
    xy.PlotLine.hide.val = not PlotLine
    x_axis = graph.x
    y_axis = graph.y
    x_axis.label.val = xname
    x_axis.log.val = xlog
    y_axis.label.val = yname
    y_axis.log.val = ylog
    #~
    #~ y_axis.min.val = -1.1
    #~ y_axis.max.val = 1.1

    #~ x_axis.min.val = 0.25
    #~ x_axis.max.val = 0.6

    xy.ErrorBarLine.width.val = '1pt'
    #~ xy_sim.ErrorBarLine.transparency.val = 50






def veusz2csv(g, filename):
    global DATA, strings, w
    # writes all the datasets from graphic object g to a csv file for further using in gnuplot
    Datasets = g.GetDatasets()
    DATA = {}
    #~ DATA = {Dataset: g.GetData(Dataset) for Dataset in Datasets}


    for Dataset in Datasets:
        DATA[Dataset] = g.GetData(Dataset)

    keys = DATA.keys()
    for i in keys:
        val = DATA[i][0]

        if not DATA[i][1] == None:
            print 'err'
            err = DATA[i][1]
            DATA[i+'_err'] = list(err)
        DATA[i] = list(val)

    import csv
    strings = []
    for i in range(len(DATA)):
        string = DATA.items()[i][1]
        string.insert(0, DATA.items()[i][0])
        strings.append(string)

    strings.sort()

    longest_length = 0

    for i in range(len(strings)):
        l = len(strings[i])
        if l > longest_length:
            longest_length = l
    for i in range(len(strings)):
        d = longest_length - len(strings[i])
        for j in range(d):
            strings[i].append(None)
    strings = array(strings)
    strings = strings.T

    with open(filename,'wb') as f:
        w = csv.writer(f, delimiter=' ')
        w.writerows(strings)
    print 'data is stored in ', filename


########Average function z#######
def calc_average(z):
 av = 0
 norm = 0

 for k in arange(0,len(z)):
  av = av+(k+1)*z[k]
  norm = norm+z[k]
 av = av/norm
 return av

########Average function z#######
def calc_2average(z):
 av2 = 0
 norm = 0

 for k in arange(0,len(z)):
  av2 = av2+(k+1)*(k+1)*z[k]
  norm = norm+z[k]
 av2 = av2/norm
 return av2


#########Fraction of stars########

def findmax(To_finding_max):
    global num
    nmax=0; l=num
    while (nmax==0) and (l-1) > 0:
	if num > 0:
         if To_finding_max[l-1] < To_finding_max[l]:
             nmax=1
             num=l
         l=l-1;
	else:
	     nmax=1
	     num=l
    return num

def findmin(To_finding_min):
    global num
    nmin=0; l=num
    while (nmin==0):
	if num > 0 and (l-1) > 0:
         if To_finding_min[l-1] > To_finding_min[l] :
            nmin=1
	    num=l
         l=l-1;
	else:
	    nmin=1
	    num=0
    return num

def lambdafind(To_finding_bp,To_finding_ne):
    global num, extremdatabp, extremdatane, lambdadictbp, lambdadictne
    netot=0; bptot=0; num=len(To_finding_ne)-1; a=0;
    extremdatabp={}; extremdatane={} ; lambdadictbp={} ; lambdadictne={}
    #Creating extremdata dictinary for branching point and ends
    num=len(To_finding_ne)-1
    j=1
    while num > 0:
	extremdatane['min'+str(j)]=findmin(To_finding_ne)
	extremdatane['max'+str(j)]=findmax(To_finding_ne)
        j=j+1
    j=1
    num=len(To_finding_bp)-1
    while num > 0:
	extremdatabp['min'+str(j)]=findmin(To_finding_bp)
	extremdatabp['max'+str(j)]=findmax(To_finding_bp)
        j=j+1
    for k in range(1,len(To_finding_ne)):
       netot = netot + To_finding_ne[k]/b.sigma/(b.f-1)
    for k in range(1,len(To_finding_bp)):
       bptot = bptot + To_finding_bp[k]/b.sigma
    mylist_bp=[]
    for key in extremdatabp.keys():
        if key[0:-1] == 'min':
           mylist_bp.append(key)
           a=len(mylist_bp)
           for i in range(1,a):
            lambdadictbp['lam'+str(i)]=calc_area_bp(To_finding_bp,extremdatabp['min'+str(i)],extremdatabp['min'+str(i+1)])
    mylist_ne=[]
    for key in extremdatane.keys():
        if key[0:-1] == 'min':
           mylist_ne.append(key)
           a=len(mylist_ne)
           for i in range(1,a):
            lambdadictne['lam'+str(i)]=calc_area_ne(To_finding_ne,extremdatane['min'+str(i)],extremdatane['min'+str(i+1)])
    return extremdatabp, extremdatane, lambdadictbp, lambdadictne

def calc_area_bp(To_finding_bp,a,c):
    lam = 0
    for k in range(a,c,-1):
        lam = lam + To_finding_bp[k]/b.sigma
    return lam

def calc_area_ne(To_finding_ne,a,c):
    lam = 0
    for k in range(a,c,-1):
        lam = lam + To_finding_ne[k]/b.sigma/(b.f-1)
    return lam

def varpar(variation):
        global my_parameter_bp, my_fraction_bp, my_parameter_ne, my_fraction_ne
        varpardictbp[variation]=lambdadictbp.values()
        varpardictne[variation]=lambdadictne.values()
        my_keys_bp=[]; my_values_bp=[]; my_keys_ne=[]; my_values_ne=[]
        for key in varpardictbp.keys():
			for length_box_value in range(0,len(varpardictbp.get(key))):
				my_keys_bp.append(key)
				my_values_bp.append(varpardictbp.get(key)[length_box_value])
        for key in varpardictne.keys():
			for length_box_value in range(0,len(varpardictne.get(key))):
				my_keys_ne.append(key)
				my_values_ne.append(varpardictne.get(key)[length_box_value])
        my_parameter_bp=np.array(my_keys_bp)
        my_fraction_bp=np.array(my_values_bp)
        my_parameter_ne=np.array(my_keys_ne)
        my_fraction_ne=np.array(my_values_ne)
        return

#########LINECOLORS AND LINETYPES ########
def my_line_color():
    my_line_color_dict={}
    my_line_color_dict={1:'black',2:'red',3:'green',4:'blue',5:'magenta',6:'#01A9DB',7:'grey',8:'#8904B1',9:'#B4045F',10:'#585858',11:'#DF3A01',12:'#DBA901',13:'#74DF00',14:'#CEF6F5'}
    #my_line_color_dict={1:'#B40404',2:'#B43104',3:'#DBA901',4:'#5FB404',5:'#01DF74',6:'#01A9DB',7:'#0040FF',8:'#8904B1',9:'#B4045F',10:'#585858'}
    line_color=my_line_color_dict[line_color_num]
    return line_color
def my_line_type():
    my_line_type_dict={}
    #my_line_type_dict={1:'solid',2:'solid',3:'solid',4:'solid',5:'solid',6:'solid',7:'solid',8:'solid',9:'solid',10:'solid'}
    #my_line_type_dict={1:'dashed',2:'dashed',3:'dashed',4:'dashed',5:'dashed',6:'dashed',7:'dashed',8:'dashed',9:'dashed',10:'dashed'}
    #my_line_type_dict={1:'dotted',2:'dotted',3:'dotted',4:'dotted',5:'dotted',6:'dotted',7:'dotted',8:'dotted',9:'dotted',10:'dotted'}
    #my_line_type_dict={1:'dash-dot',2:'dash-dot',3:'dash-dot',4:'dash-dot',5:'dash-dot',6:'dash-dot',7:'dash-dot',8:'dash-dot',9:'dash-dot',10:'dash-dot'}
    my_line_type_dict={1:'solid',2:'dashed',3:'dotted',4:'dash-dot',5:'dash-dot-dot',6:'dotted-fine',7:'dashed-fine',8:'dash-dot-fine',9:'dot1',10:'dot2'}
    line_type=my_line_type_dict[line_type_num]
    return line_type

#########FRACTION OF EXTENDED STARS ########
def extendedStarFind(To_finding_bp,To_finding_ne):
    global num, firstdatabp, firstdatane, lambdaonebp, lambdaonene
    netot=0; bptot=0; num=len(To_finding_ne)-1; a=0;
    firstdatabp={}; firstdatane={} ; lambdaonebp={} ; lambdaonene={}
    #Creating e:xtremdata dictinary for branching point and ends
    num=len(To_finding_ne)-1
    j=1
    while j < 3:
	firstdatane['min'+str(j)]=findmin(To_finding_ne)
	firstdatane['max'+str(j)]=findmax(To_finding_ne)
        j=j+1
    j=1
    num=len(To_finding_bp)-1
    while j < 3:
	firstdatabp['min'+str(j)]=findmin(To_finding_bp)
	firstdatabp['max'+str(j)]=findmax(To_finding_bp)
        j=j+1
    for k in range(1,len(To_finding_ne)):
       netot = netot + To_finding_ne[k]/b.sigma/(b.f-1)
    for k in range(1,len(To_finding_bp)):
       bptot = bptot + To_finding_bp[k]/b.sigma
    if firstdatabp.get('min2') == 0:
       lambdaonebp['lam1']=0
    else:
       lambdaonebp['lam1']=calc_area_bp(To_finding_bp,firstdatabp['min1'],firstdatabp['min2'])
    if firstdatane.get('min2') == 0:
       lambdaonene['lam1']=0
    else:
       lambdaonene['lam1']=calc_area_ne(To_finding_ne,firstdatane['min1'],firstdatane['min2'])
    return firstdatabp, firstdatane, lambdaonebp, lambdaonene

def varparfirst(variation):
        global my_oneparameter_bp, my_onefraction_bp, my_oneparameter_ne, my_onefraction_ne, my_twoparameter_bp, my_twoparameter_ne
        varparfirstbp[variation]=lambdaonebp.values()
        varparfirstne[variation]=lambdaonene.values()
        my_keys_bp=[]; my_values_bp=[]; my_keys_ne=[]; my_values_ne=[]; my_keys2_bp=[]; my_keys2_ne=[]
        for key in varparfirstbp.keys():
			for length_box_value in range(0,len(varparfirstbp.get(key))):
				my_keys_bp.append(key)
				my_keys2_bp.append(b.chi)
				my_values_bp.append(varparfirstbp.get(key)[length_box_value])
        for key in varparfirstne.keys():
			for length_box_value in range(0,len(varparfirstne.get(key))):
				my_keys_ne.append(key)
				my_keys2_ne.append(b.chi)#####
				my_values_ne.append(varparfirstne.get(key)[length_box_value])
        my_oneparameter_bp=np.array(my_keys_bp)
        my_twoparameter_bp=np.array(my_keys2_bp)
        my_onefraction_bp=np.array(my_values_bp)
        my_oneparameter_ne=np.array(my_keys_ne)
        my_twoparameter_ne=np.array(my_keys2_ne)
        my_onefraction_ne=np.array(my_values_ne)
        return


#########CALCULATION PROPAGATOR GTl AND GF (LINE CASE)########
def calc_gtl_gfl(w):
    global gtl, gfl
    #init
    gtl = np.zeros((b.chain_length,b.chain_length+1))
    gfl = np.zeros((b.chain_length,b.chain_length+1))
    #condition
    gtl[0][0] = w[0]
    for j in range(0,b.chain_length+1):
        gfl[0][j] = w[j]
    #start_recurrence
    for k in range(1,b.chain_length):
        gtl[k][0] = num_lambda*w[0]*(4.0*gtl[k-1][0]+gtl[k-1][1])
        gfl[k][0] = num_lambda*w[0]*(4.0*gfl[k-1][0]+gfl[k-1][1])
   	for j in range(1,k+1):
            gtl[k][j] = num_lambda*w[j]*(gtl[k-1][j-1]+4.0*gtl[k-1][j]+gtl[k-1][j+1])
        for j in range(1,b.chain_length-k+1):
            gfl[k][j] = num_lambda*w[j]*(gfl[k-1][j-1]+4.0*gfl[k-1][j]+gfl[k-1][j+1])
    return gtl, gfl

#########CALCULATION PHI AND PROBABILYTI OF TERMINAL GROUP (LINE CASE)########
def calc_phiL_zeL(w):
    global phiL, zeL
    Zn = gfl[b.chain_length-1][0]
    phiL = np.zeros((b.chain_length))
    zeL = np.zeros((b.chain_length))
    for k in range(0,b.chain_length-1):
        for j in range(0,b.chain_length-1):
            phiL[k] = phiL[k] + gtl[j][k]*gfl[b.chain_length-j-1][k]/w[k]
    phiL = phiL/Zn
    phiL = phiL/b.chain_length
    zeL = gtl[b.chain_length-1]/Zn
    return phiL, zeL

#########CALCULATION PROPAGATOR GT AND GF (BRUSH CASE)########
def calc_gtb_gfb(w):
    global gtb, gfb
    #init
    gtb = np.zeros((arm_length+2,arm_length+2))
    gfb = np.zeros((arm_length+2,2*arm_length+2))
    #condition
    gtb[0][0] = w[0]
    for j in range(0,(2*arm_length+1)+1):
        gfb[0][j] = w[j]
    #start_recurrence
    for k in range(1,arm_length+1):
        gtb[k][0] = num_lambda*w[0]*(4.0*gtb[k-1][0]+gtb[k-1][1])
        gfb[k][0] = num_lambda*w[0]*(4.0*gfb[k-1][0]+gfb[k-1][1])
        for j in range(1,k+1):
            gtb[k][j] = num_lambda*w[j]*(gtb[k-1][j-1]+4.0*gtb[k-1][j]+gtb[k-1][j+1])
        for j in range(1,(2*arm_length+1)-k-1):
            gfb[k][j] = num_lambda*w[j]*(gfb[k-1][j-1]+4.0*gfb[k-1][j]+gfb[k-1][j+1])
    return gtb, gfb

#########CALCULATION PROPAGATOR G2 (BRUSH CASE)########
def calc_g2(w):
    global g2
    #init
    g2 = np.zeros((arm_length+2,arm_length+2,2*arm_length+5))
    #condition
    for j in range(0,arm_length+1):
        g2[0][j][j] = w[j]
    #start_recurrence
    for k in range(1,arm_length+1):
        for z1 in range(0,arm_length+1):
            g2[k][z1][0] = num_lambda*w[0]*(4.0*g2[k-1][z1][0]+g2[k-1][z1][1])
            for j in range(1,k+z1+1):
                g2[k][z1][j] = num_lambda*w[j]*(g2[k-1][z1][j-1]+4.0*g2[k-1][z1][j]+g2[k-1][z1][j+1])
    return g2

#########CALCULATION ZBP AND ZE (BRUSH CASE)########
def calc_zbp_ze(w):
    global zbp_calc, ze_calc
    Znb = 0; Zne = 0
    zbp_calc = np.zeros((arm_length+2))
    ze_calc = np.zeros((2*arm_length+1))
    for k in range(0,arm_length):
        zbp_calc[k] = gtb[arm_length-1][k]*(gfb[arm_length][k]/w[k])**(f1-1)
        Znb = Znb + zbp_calc[k]
    zbp_calc = zbp_calc/Znb
    for k in range(0,2*arm_length+1):
        for z1 in range(0,arm_length):
            ze_calc[k] = ze_calc[k] + gtb[arm_length-1][z1]*(gfb[arm_length][z1]/w[z1])**(f1-2)*(g2[arm_length][z1][k]/w[z1])
        Zne = Zne + ze_calc[k]
    ze_calc = ze_calc/Zne
    return zbp_calc, ze_calc

def potention_from_formula():
    global M
    M = (math.pi*n_formula/2)/arccos(sqrt((f_formula-1.)/f_formula))
    potention_formula = np.zeros((b.chain_length+1))
    for i in range(0,b.chain_length+1):
        if i < (b.h2n*2*n_formula-1):
            #potention_formula[i] = -3*log(cos(b.h2n*math.pi/2))-((3*math.pi**2/8)*((b.h2n)**2-(b.h2n*2*n_formula/M)**2))-(-3*log(cos((math.pi/2)*(i/200)))-((3*math.pi**2/8)*((i/200)**2-(i/M)**2)))
            potention_formula[i] = -3.*log(cos(b.h2n*math.pi/2.))-((3.*math.pi**2/8.)*((b.h2n)**2-(b.h2n*2.*n_formula/M)**2))-(-3.*log(cos((math.pi/2.)*(i/(2.*n_formula))))-((3.*math.pi**2/8.)*((i/(2.*n_formula))**2-(i/M)**2)))
        else:
            potention_formula[i] = 0
    return potention_formula

def change_potention_from_sfbox(pot_in):
    global M
    pot_up = np.zeros((2*n_formula+5))
    pot_mod = np.zeros((2*n_formula+5))
    pot_out = np.zeros((b.chain_length+5))
    M = (math.pi*n_formula/2)/arccos(sqrt((f_formula-1.)/f_formula))
    pot_up = pot_in[5] - pot_in
    for i in range(0,(2*n_formula)):
        pot_mod[i] = pot_up[i] -((3.*math.pi**2/8.)*((i/(2.*n_formula))**2-(i/M)**2))
    for i in range(0,b.chain_length+5):
        if i < (b.h2n*2*n_formula-1):
            pot_out[i] = pot_mod[b.h2n*2*n_formula-1]-pot_mod[i]
        else:
            pot_out[i] = 0
    pot_out = array(pot_out)
    return pot_out

#~ for p in range(len(To_finding_h),0,-1):
        #~ if  To_finding_h[p-1]<1.0e-8:
         #~ h=p
        #~ p=p-1;
    #~ h1=min(b.n, extremdata['min2'])
    #~ h2=h-h1
