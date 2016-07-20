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


#Average function z
def z_av(z):
 z_av = 0
 z2_av = 0
 norm = 0

 for k in arange(0,len(z)):
  z_av = z_av+k*z[k]
  z2_av = z2_av+k*k*z[k]
  norm = norm+z[k]
 z_av = z_av/norm
 return z


#########Fraction of extended stars########   

def max(To_finding_max):
    global num
    nmax=0 
    for l in range(num,0,-1):
            while (nmax==0):
                if To_finding_max[l-2] < To_finding_max[l-1]:
                   nmax=1
                   num=l-1
                l=l-1;
    return num

def min(To_finding_min):
    global num
    nmin=0
    for l in range(num,0,-1):
            while (nmin==0):
                if To_finding_min[l-2] > To_finding_min[l-1] :
                  nmin=1
                  num=l-1
                l=l-1;
    return num

def lambdafind(To_finding_bp,To_finding_ne):
    global num, extremdatabp, extremdatane, lambdadictbp, lambdadictne
    netot=0; bptot =0; num=len(To_finding_ne); a = 0;
    extremdatabp={}; extremdatane={} ; lambdadictbp={} ; lambdadictne={}
    j=1
    #Creating extremdata dictinary for branching point and ends
    while num>0:
     if (min(To_finding_ne))>0 or (max(To_finding_ne))>0:
      extremdatane['min'+str(j)]=min(To_finding_ne)
      extremdatane['max'+str(j)]=max(To_finding_ne)
     if (min(To_finding_bp))>0 or (max(To_finding_bp))>0:
      extremdatabp['min'+str(j)]=min(To_finding_bp)
      extremdatabp['max'+str(j)]=max(To_finding_bp)
     j=j+1
    extremdatabp['min'+str(j-1)]=0
    extremdatane['min'+str(j-1)]=0
    for k in range(1,len(To_finding_ne)):
       netot = netot + To_finding_ne[k]/b.sigma/(b.f-1)
    for k in range(1,len(To_finding_bp)):
       bptot = bptot + To_finding_bp[k]/b.sigma
    mylist_ne=[]
    for key in extremdatane.keys():
        if key[0:-1] == 'min':
           mylist_ne.append(key)
           a=len(mylist_ne)
           for i in range(1,a):
            lambdadictne['lam'+str(i)]=calc_area_ne(To_finding_ne,extremdatane['min'+str(i)],extremdatane['min'+str(i+1)])
    mylist_bp=[]
    for key in extremdatabp.keys():
        if key[0:-1] == 'min':
           mylist_bp.append(key)
           a=len(mylist_bp)
           for i in range(1,a):
            lambdadictbp['lam'+str(i)]=calc_area_bp(To_finding_bp,extremdatabp['min'+str(i)],extremdatabp['min'+str(i+1)])
    return lambdadictbp, lambdadictne
    
   
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
    
    #~ for p in range(len(To_finding_h),0,-1):
        #~ if  To_finding_h[p-1]<1.0e-8:
         #~ h=p
        #~ p=p-1;
    #~ h1=min(b.n, extremdata['min2'])
    #~ h2=h-h1
           


#~ def lambdafind(To_finding_tp,To_finding_h):
    #~ #lambdafile = open('data/'+self.fnameout)
    #~ #phi_tolambda = array(self.datadict['mol : brush : phi'])
    #~ #ne_tolambda = array(self.datadict['mol : brush : phi-polyn'])
 #~ 
    #~ nmax=0 ; nmin=0 ; h=0;  lam=0;  netot=0;
    #~ 
    #~ for l in range(len(To_finding_tp),0,-1):
            #~ while (nmax==0):
                #~ if To_finding_tp[l-2] < To_finding_tp[l-1] :
                   #~ nmax=1
                   #~ lmax=l-1
                #~ l=l-1;
    #~ for k in range(lmax,0,-1):
            #~ while (nmin==0):
                #~ if To_finding_tp[k-2] > To_finding_tp[k-1] :
                  #~ nmin=1
                  #~ kmin=k-1
                #~ k=k-1;
    #~ for p in range(len(To_finding_h),0,-1):
        #~ if  To_finding_h[p-1]<1.0e-8:
         #~ h=p
        #~ p=p-1;
    #~ h1=min(b.n, kmin)
    #~ h2=h-h1
    #~ for k in range(kmin,len(To_finding_tp)):
        #~ lam = lam + To_finding_tp[k]/b.sigma/(b.f-1)
    #~ 
    #~ for k in range(1,len(To_finding_tp)):
        #~ netot = netot + To_finding_tp[k]/b.sigma/(b.f-1)
    #~ return lmax, kmin, h, h1, h2, lam/netot
