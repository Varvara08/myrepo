# -*- coding: utf-8 -*-
# This is program for creating .in file for sfbox and calculating it.
# The brush is consist of linear chain --- simple linear brush (f=2).
# Ploting by veusz.embeded
#
#!/usr/bin/python
from scipy import *
from scipy.optimize import fsolve
from pylab import *
import numpy as np
import os

from copy import deepcopy
from io import StringIO
import veusz.embed as veusz
execfile("functions.py")

import sys

class brush():
    isPlusOne = False; # Rbox = Rstar/A + isPlusOne

    val = -1; # valency of monomers
    #Rstar = 10; # Radius at which the end monomers are pinned
    #~ Rbox = 20;
    n = 100 # number of segments per arm
    f = 3 # number of arms
    sigma = 0.1  # density of grafting brush
    line_sigma = 1e-6  # density of grafting chain
    chi = 0.0
    pK = 6.0
    epsilon = 80
    chain_length = 100
    cl=1e-1 # phibulk(chloride)
    pH=7 # (unused)
    iguess_in=''; fname_prefix='';

    solved = False
    # A constants
    #~ distance = 3.1044e-10 # [m] The space unity in sfbox model
    distance = 3.6e-10 # [m] The space unity in sfbox model
    Navogadro = 6.022e23 # [1/mol] Avogadro number
    lB = 0.7e-9
    kT=4.14e-21; # [J] 1kT = 4.14e-21 J at 300K
    to_mols = 1/Navogadro/ distance**3/1000 # [mol/l]

    def __init__(self):
        self.rD = int( (8*pi*self.lB*self.Navogadro*self.cl*self.to_mols / 0.001)**-0.5 / self.distance)
        self.Zbox =  int(self.n*2+2*self.rD)
        self.fname="brush_d+l_n"+str(self.n)+"_f"+str(self.f)+"_chain_length"+str(self.chain_length)+"_val"+str(self.val)+'_cl{:.2e}'.format(self.cl)+'_chi'+str(self.chi)+"_pK"+str(self.pK)+"_pH"+str(self.pH)+'_sig'+str(self.sigma)
        self.fnamein = self.fname+".in"
        self.fnameout = self.fname+".out"

    def inputfilegen(self):
        self.__init__()
        n = self.n
        chi = self.chi; pK=self.pK; cl = self.cl
        Zbox = self.Zbox

        output=[]
        output.append('')

        output.append('mon : na : valence : 1')
        output.append('mon : cl : valence : -1')

        output.append('mol : na : composition : na')
        output.append('mol : na : freedom : neutralizer')

        output.append('mol : cl : composition : cl')
        output.append('mol : cl : freedom : free')
        output.append('mol : cl : phibulk : '+str(self.cl))

        output.append('mol : solvent : composition : water')
        output.append('mol : solvent : freedom : solvent')

        segments = ['poly','poly1','polyE','polyB']
        freedom  = ['free','pinned', 'free','free']

        #~ freedom  = ['free','pinned','free','free','pinned']

        for i in range(len(segments)):
            output.append('mon : '+segments[i]+' : freedom : '+freedom[i])

        output.append('mon : water : freedom : free')

        composition = '(poly1)'+'(poly)'+str(self.n-2)+'(polyB)'+(self.f-2)*('[(poly)'+str(self.n-1)+'(polyE)]') +'(poly)'+str(self.n-1)+'(polyE)'
        self.composition = composition
        output.append('mol : brush : composition : '+composition)

	output.append('//////////////// LINE ///////////////')

	line_segments = ['Lpoly','Lpoly1','LpolyE']
        line_freedom  = ['free','pinned', 'free']

	for i in range(len(line_segments)):
            output.append('mon : '+line_segments[i]+' : freedom : '+line_freedom[i])

        composition = '(poly1)'+'(poly)'+str(self.chain_length-2)+'(polyE)'
        self.composition = composition
        output.append('mol : line : composition : '+composition)


        output.append('mol : line : freedom : restricted')
        output.append('mol : line : theta : '+str(self.chain_length*self.line_sigma))

        output.append('mon : Lpoly1 : pinned_range : 1;1')
	output.append('//////////////// LINE END///////////////')

	if self.pK != 0:
            import fractions
            alpha = abs(self.val)
            alpha_num = fractions.Fraction(alpha).numerator
            alpha_den = fractions.Fraction(alpha).denominator
            output.append('//////////////// Ionization ///////////////')
            output.append('mon : water : state1 : H3O')
            output.append('mon : water : state2 : H2O')
            output.append('mon : water : state3 : OH')
            output.append('state : H3O : valence : 1')
            output.append('state : H3O : alphabulk : '+str(10**(-self.pH)))
            output.append('state : H2O : valence : 0')
            output.append('state : OH : valence : -1')
            output.append('reaction : auto : equation : 2(H2O)=1(H3O)1(OH)')
            output.append('reaction : auto : pK : 17.5')
            # ?? make segment+HA a new string variable
            for segment in segments:
                output.append('mon : '+segment+' : state1 : '+segment+'HA')
                output.append('mon : '+segment+' : state2 : '+segment+'A')
                output.append('state : '+segment+'HA : valence : 0')
                output.append('state : '+segment+'A : valence : '+str(self.val))
                output.append('reaction : weak_'+segment+' : equation : '+str(alpha_den)+'('+segment+'HA)'+str(alpha_num)+'(H2O)='+str(alpha_den)+'('+segment+'A)'+str(alpha_num)+'(H3O)')
                output.append('reaction : weak_'+segment+' : pK : '+str(pK))
            for line_segment in line_segments:
                output.append('mon : '+line_segment+' : state1 : '+line_segment+'HA')
                output.append('mon : '+line_segment+' : state2 : '+line_segment+'A')
                output.append('state : '+line_segment+'HA : valence : 0')
                output.append('state : '+line_segment+'A : valence : '+str(self.val))
                output.append('reaction : weak_'+line_segment+' : equation : '+str(alpha_den)+'('+line_segment+'HA)'+str(alpha_num)+'(H2O)='+str(alpha_den)+'('+line_segment+'A)'+str(alpha_num)+'(H3O)')
                output.append('reaction : weak_'+line_segment+' : pK : '+str(pK))
        else:
            for segment in segments:
                output.append('mon : '+segment+' : valence : '+str(self.val))
            for line_segment in line_segments:
                output.append('mon : '+line_segment+' : valence : '+str(self.val))


        output.append('/////////////////// Chi ///////////////////')
        for segment in segments:
            output.append('mon : '+segment+' : chi - water : '+str(chi))
            output.append('mon : '+segment+' : chi - na : '+str(chi))
            output.append('mon : '+segment+' : chi - cl : '+str(chi))
        for line_segment in line_segments:
            output.append('mon : '+line_segment+' : chi - water : '+str(chi))
            output.append('mon : '+line_segment+' : chi - na : '+str(chi))
            output.append('mon : '+line_segment+' : chi - cl : '+str(chi))


        output.append('////////////////////////////////////////////')

        output.append('mol : brush : freedom : restricted')
        output.append('mol : brush : theta : '+str(self.n*self.f*self.sigma))

        output.append('lat : mylat : geometry : flat')
        output.append('lat : mylat : n_layers : '+str(Zbox))
        output.append('lat : mylat : lambda : '+str(num_lambda))
        output.append('lat : mylat : lowerbound : surface')
        output.append('mon : S : freedom : frozen')
        output.append('mon : S : frozen_range : lowerbound')
        output.append('/////////////////////////////////////////////')

        output.append('output : filename.out : type : ana')
        output.append('output : filename.out : write_bounds : false')

        output.append('/////////////////////////////////////////////')
        output.append('sys : noname : overflow_protection : true')
        output.append('newton : isac : iterationlimit : 100000')
        output.append('newton : isac : tolerance : 1e-7')
        #~ output.append('//newton : isac : method : LBFGS')
        #~ output.append('//newton : isac : m : 8')
        self.iguess_out=self.fname+'.iguess'
        output.append('newton : isac : initial_guess_output_file : '+self.iguess_out)
        if self.iguess_in != '':
            if os.path.exists(mypath+'data/'+self.iguess_in):
                output.append('newton : isac : initial_guess : file')
                output.append('newton : isac : initial_guess_input_file : '+self.iguess_in)
            else:
                self.iguess_in = ''
        else:
            #~ output.append('//newton : isac : initial_guess_output_file : X')
            output.append('//newton : isac : initial_guess : file')
            output.append('//newton : isac : initial_guess_input_file : X')

        self.iguess_in = self.iguess_out
        output.append('/////////////////////////////////////////////')

        output.append('mon : poly1 : pinned_range : 1;1')

        inputfile = open(mypath+'data/'+self.fnamein,'w')
        #~ print (self.fnamein)
        for i in range(len(output)):
            inputfile.writelines(output[i]+'\n');
        inputfile.close()
        print self.fname
        return output

    def sfbox(self):
        global sfbox_out
        self.solved = False
        self.__init__()
        if os.path.exists(mypath+'data/'+self.fnameout):
            #~ print ("task is already solved")
            self.solved = True
            self.loadData()
            #~ self.Nna = float(self.getValue('mol : na : theta excess : '))
            #~ self.Ncl = float(self.getValue('mol : cl : theta excess : '))

        else:
            #~ os.system('mv '+fnametmp+" "+fnamein)
            #~ sfbox_out = os.popen("sfbox "+fnamein+' | tee /dev/pts/0').read()
            infile = self.inputfilegen()
            os.chdir(mypath+'data')

            sfbox_out = os.popen("sfbox "+self.fnamein).read()
            os.chdir('..')
            sfbox_out = sfbox_out.split("\n")

            print (sfbox_out[-3])
            if sfbox_out[-3] == 'Problem solved.':
                self.solved = True
                self.loadData()
                print (self.fnamein)
            else:
                os.remove(mypath+'data/'+self.fnameout)
                self.iguess_in = ''
                print ('Attempt number two')
                infile = self.inputfilegen()
                os.chdir(mypath+'data')
                sfbox_out = os.popen("sfbox "+self.fnamein).read()
                os.chdir('..')
                sfbox_out = sfbox_out.split("\n")
                print (sfbox_out[-3])
                if sfbox_out[-3] == 'Problem solved.':
                    self.solved = True
                    self.loadData()

                else:
                    self.solved = False
                    os.remove(mypath+'data/'+self.fnameout)
                print (self.fnamein)

            print
        return self.solved

    def loadData(self):
        global datadict, datastring, i, data, value
        datadict = {}
        datafname = mypath+'data/'+self.fnameout
        datafile = open(datafname, 'rb');
        data = datafile.read(); datafile.close()
        data = data.split('\n')
        data = data[:-3]

        keys = []
        values = []
        i = 0
        while i < len (data)-1:
            #~ print data[i].split(' : ')[-1]
            if(data[i]=="system delimiter" or len(data[i])==0):
                i = i+1;
                continue;
            elif data[i].split(' : ')[-1] == 'profile' or data[i].split(' : ')[-1] == 'vector':

                #~ key = data[i]
                key = data[i][:data[i].rindex(':')-1]
                i = i+1
                value = []
                while True:
                    try:
                        value.append(float(data[i]))
                        i = i+1
                    except (ValueError, IndexError): break
                datadict[key] = value

            else:
                key = data[i][:data[i].rindex(':')-1]
                try: value = float(data[i].split(' : ')[-1])
                except ValueError: value = data[i].split(' : ')[-1]
                datadict[key] = value
                i = i+1

        self.datadict = datadict
        return datadict

    def mylambda(self):
        global lambdadictbp, lambdadictne, ze, zbp
        zbp = array(self.datadict['mol : brush : phi-polyB'])
        ze = array(self.datadict['mol : brush : phi-polyE'])
        lambdafind(zbp,ze)
        extendedStarFind(zbp,ze)
        varpar(vary_parameter)
        varparfirst(vary_parameter)
        return lambdadictbp, lambdadictne, ze, zbp

    def chainInPotention(self):
        if (len(array(self.datadict['mol : brush : phi']))) < (chain_length+5):
          phi_list = self.datadict['mol : brush : phi']
          for k in range(0,chain_length):
            phi_list.append(0.0)
        else:
          phi_list = self.datadict['mol : brush : phi']
        phi = array(phi_list)
        potention = -log(1-phi)
        w = exp(-potention)
        print 'Calculation gtl and gfl...'
        calc_gtl_gfl(w)
        print 'Complete.'
        calc_phiL_zeL(w)

    def brushInPotention(self):
        global w
        phi= array(self.datadict['mol : brush : phi'])
        potention_from_formula()
        potention = -log(1-phi)
        w = exp(-potention)
        print 'Calculation gtb and gfb...'
        calc_gtb_gfb(w)
        print 'Complete.'
        print 'Calculation g2...'
        calc_g2(w)
        print 'Complete.'
        calc_zbp_ze(w)

    def plotProfile(self, profile_name = 'sys : noname : potential'):
        global g, xy, y_axis
        line_color=my_line_color()
        #line_type=my_line_type()
        point_color_bp='#'+3*(hex(int((1-my_onefraction_bp[0])*255))[-2:])
        point_color_ne='#'+3*(hex(int((1-my_onefraction_ne[0])*255))[-2:])
        markerSize = '1pt'
        dataprefix = self.fname+'_'
        #~ dataprefix = ''
        try:
            type (g) == veusz.Embedded
            page = g.Root.page1
            graph0 = g.Root.page1.grid1.graph1
            graph1 = g.Root.page1.grid1.graph2
            graph2 = g.Root.page1.grid1.graph3
            graph3 = g.Root.page1.grid1.graph4
            graph4 = g.Root.page1.grid1.graph5
            graph5 = g.Root.page1.grid1.graph6
        except (NameError, AttributeError):
            g = veusz.Embedded(self.fname)
            g.EnableToolbar()

            page = g.Root.Add('page')
            page.width.val = '20cm'
            page.height.val = '20cm'
            g.ResizeWindow(800,900)
            grid = page.Add('grid')
            grid.bottomMargin.val ='0cm'
            grid.leftMargin.val = '0cm'
            grid.rightMargin.val = '0cm'
            grid.topMargin.val = '0cm'
            grid.rows.val = 4
            grid.columns.val = 2
            graph0 = grid.Add('graph')
            graph1 = grid.Add('graph')
            graph2 = grid.Add('graph')
            graph3 = grid.Add('graph')
            graph4 = grid.Add('graph')
            graph5 = grid.Add('graph')

        #~ g.To('/page1/grid1/graph1')

        self.loadData()
        z = arange(0,self.Zbox)
        z = z[1:]
        x_dataname = dataprefix+'z'
        g.SetData(x_dataname, z)

        ########################
        graph = graph0
        ## Plot polymer dencity #
        what_toPlot = 'phi'
        profile_toPlot = array(self.datadict['mol : brush : phi'])
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = line_color
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = line_color

        y_axis=graph.y
        x_axis=graph.x

        #y_axis.log.val = True
        y_axis.min.val = 1e-4
        x_axis.max.val = 200
        y_axis.label.val = '\italic{\\phi}'
        x_axis.label.val = '\italic{z}'

        ########################
        graph = graph1
        ## Plot counterions #
        what_toPlot = 'phi_Na'
        profile_toPlot = array(self.datadict['mol : na : phi'])
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = line_color
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = line_color

        ########################
        graph = graph1
        ## Plot coions #
        what_toPlot = 'phi_Cl'
        profile_toPlot = array(self.datadict['mol : cl : phi'])
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = line_color
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = line_color

        y_axis=graph.y
        x_axis=graph.x

        y_axis.log.val = True
        #y_axis.min.val = 1e-9
        #y_axis.max.val = 1e-1
        #x_axis.max.val = 250
        y_axis.label.val = '\italic{\\phi_{Na}, \\phi_{Cl}}'
        x_axis.label.val = '\italic{z}'

        ########################
        graph = graph3
        # Plot the distribution of the ends
        what_toPlot = 'ends'
        profile_toPlot = array(self.datadict['mol : brush : phi-polyE'])
        profile_toPlot=profile_toPlot/b.sigma/(b.f-1)
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = line_color
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = line_color

        y_axis=graph.y
        x_axis=graph.x

        #y_axis.log.val = True
        #y_axis.min.val = 1e-8
        x_axis.max.val = 200
        y_axis.label.val = '\italic{n_{ends}}'
        x_axis.label.val = '\italic{z}'

        ########################
        graph = graph2
        # Plot the distribution of the branching points
        what_toPlot = 'bp'
        profile_toPlot = array(self.datadict['mol : brush : phi-polyB'])
        profile_toPlot=profile_toPlot/b.sigma
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = line_color
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = line_color

        y_axis=graph.y
        x_axis=graph.x

        #y_axis.log.val = True
        #y_axis.min.val = 1e-8
        x_axis.max.val = 110
        y_axis.label.val = '\italic{n_{bp}}'
        x_axis.label.val = '\italic{z}'

        ########################
        graph = graph4
        # Plot chain dencity
        what_toPlot = 'line_phi_sfbox'
        profile_toPlot = array(self.datadict['mol : line : phi'])
        profile_toPlot=profile_toPlot/(self.chain_length*self.line_sigma)
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = line_color
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = line_color

        y_axis=graph.y
        x_axis=graph.x

        #y_axis.log.val = True
        #y_axis.min.val = 1e-4
        #x_axis.max.val = 200
        y_axis.label.val = '\italic{\\phi_l}'
        x_axis.label.val = '\italic{z}'


        ########################
        graph = graph5
        # Plot the distribution of the chain ends (SFBox)
        what_toPlot = 'line_ends_sfbox'
        profile_toPlot = array(self.datadict['mol : line : phi-polyE'])
        profile_toPlot=profile_toPlot/b.line_sigma
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = line_color
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = line_color

        y_axis=graph.y
        x_axis=graph.x

        #y_axis.log.val = True
        #y_axis.min.val = 1e-8
        #x_axis.max.val = 200
        y_axis.label.val = '\italic{nl_{ends}}'
        x_axis.label.val = '\italic{z}'


    def plotDiagram(self, profile_name = 'sys : noname : potential'):
        global g, xy, y_axis
        #line_color=my_line_color()
        #line_type=my_line_type()
        point_color_bp='#'+3*(hex(int((1-my_onefraction_bp[0])*255))[-2:])
        point_color_ne='#'+3*(hex(int((1-my_onefraction_ne[0])*255))[-2:])
        markerSize = '1pt'
        dataprefix = self.fname+'_'
        #~ dataprefix = ''
        try:
            type (g) == veusz.Embedded
            page = g.Root.page1
            graph0 = g.Root.page1.grid1.graph1
            graph1 = g.Root.page1.grid1.graph2
        except (NameError, AttributeError):
            g = veusz.Embedded(self.fname)
            g.EnableToolbar()

            page = g.Root.Add('page')
            page.width.val = '20cm'
            page.height.val = '10cm'
            g.ResizeWindow(900,600)
            grid = page.Add('grid')
            grid.bottomMargin.val ='0cm'
            grid.leftMargin.val = '0cm'
            grid.rightMargin.val = '0cm'
            grid.topMargin.val = '0cm'
            grid.rows.val = 1
            grid.columns.val = 2
            graph0 = grid.Add('graph')
            graph1 = grid.Add('graph')

        ########################
        # Plot the extended fractions calc on branching point ###
        x_dataname = dataprefix+'pH_bp'
        g.SetData(x_dataname, my_twoparameter_bp)
        # Plot the fraction of extended stars
        graph = graph0
        ########################
        markerSize = '2pt'
        # Plot lambda#
        what_toPlot = 'lambdaone_bp'
        profile_toPlot = my_oneparameter_bp # my_oneparameter_bp --- it's already array
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'square'

        xy.MarkerFill.color.val = point_color_bp
        xy.MarkerLine.color.val = '#eeeeee'
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '2pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.hide.val = not False # False <---> don't plot lines

        y_axis=graph.y
        x_axis=graph.x

#        y_axis.log.val = True
        y_axis.min.val = 1e-4
        y_axis.max.val = 1e-1
        x_axis.min.val = 2.0
        x_axis.max.val = 15.0
        y_axis.label.val = '\italic{\\sigma}'
        x_axis.label.val = '\italic{pH}'

        ########################
        # Plot the extended fractions calc on terminal group
        x_dataname = dataprefix+'pH_ne'
        g.SetData(x_dataname, my_twoparameter_ne)
        # Plot the fraction of extended stars
        graph = graph1
        ########################
        markerSize = '2pt'
        # Plot lambda#
        what_toPlot = 'lambdaone_ne'
        profile_toPlot = my_oneparameter_ne # my_oneparameter_ne --- it's already array
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'square'

        xy.MarkerFill.color.val = point_color_ne
        xy.MarkerLine.color.val = '#eeeeee'
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '2pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.hide.val = not False # False <---> don't plot lines

        y_axis=graph.y
        x_axis=graph.x

#        y_axis.log.val = True
        y_axis.min.val = 1e-4
        y_axis.max.val = 1e-1
        x_axis.min.val = 2.0
        x_axis.max.val = 15.0
        y_axis.label.val = '\italic{\\sigma}'
        x_axis.label.val = '\italic{pH}'

    def getFreeEnergy(self):
        F = self.datadict['sys : noname : Helmholtz energy(po)']
        #~ F = self.datadict['sys : noname : free energy']
        #~ mu1 = self.datadict['mol : na : sum n FH-MU']
        #~ mu2 = self.datadict['mol : cl : sum n FH-MU']
        #~ mu3 = self.datadict['mol : solvent : sum n FH-MU']
        #~ F = F - mu1 - mu2 - mu3

        #~ sites = self.datadict['lat : mylat : lattice sites']
        #~ nlayers = self.datadict['lat : mylat : n_layers']
        #~ V = self.datadict['lat : mylat : total number of lattice sites']
        #~ for r in arange(self.Rstar, self.Rbox+0):
            #~ mu = self.datadict['mol : medium'+str(r)+' : FH-MU']
            #~ theta = self.datadict['mol : medium'+str(r)+' : theta']
            #~ entropy = self.datadict['mol : medium'+str(r)+' : entropy']
            #~ entropy_sites = entropy/sites[r]

            #~ mu = (mu + entropy/theta)*theta
            #~ entropy = entropy*V
            #~ print 'F  = ',F;
            #~ print 'mu = ',mu;
            #~ print 's  = ',entropy;
            #~ F = F - mu

            #~ print 'entropy = ', entropy/sites[r]
        #~ print self.fname
        #~ print self.Rstar
        #~ print self.Rbox
        #~ print nlayers
        #~ F= F - entropy_sites*V
        return F
