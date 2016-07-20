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
    n = 100; # number of segments per arm
    f = 3; # number of arms
    sigma = 0.1;  # density of grafting
    chi=0.0;
    pK=6.0;
    epsilon=80;
    cl=1e-1; # phibulk(chloride)
    pH=7; # (unused)
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
        #~ self.chi = chi
        self.Zbox =  int(self.n*self.f*1.5)
        self.rD = int( (8*pi*self.lB*self.Navogadro*self.cl*self.to_mols / 0.001)**-0.5 /self.distance)
        self.fname='brush_d_Z'+str(self.Zbox)+"_n"+str(self.n)+"_f"+str(self.f)+"_val"+str(self.val)+'_cl{:.2e}'.format(self.cl)+'_chi'+str(self.chi)+"_pK"+str(self.pK)+'_sig'+str(self.sigma)
        self.fnamein = self.fname+".in"
        self.fnameout = self.fname+".out"
#        self.fnamedat = self.fname+".dat" # for two population of brush

    def inputfilegen(self):
        self.__init__()
        n = self.n
        chi = self.chi; pK=self.pK; cl = self.cl
        Zbox = self.Zbox # probably here would put rd, but I'm not sure

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
            output.append('state : H3O : alphabulk : 1e-7')
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

        else:
            for segment in segments:
                output.append('mon : '+segment+' : valence : '+str(self.val))


        output.append('/////////////////// Chi ///////////////////')
        for segment in segments:
            output.append('mon : '+segment+' : chi - water : '+str(chi))
            output.append('mon : '+segment+' : chi - na : '+str(chi))
            output.append('mon : '+segment+' : chi - cl : '+str(chi))


        output.append('////////////////////////////////////////////')

        output.append('// n=30, f=5')


        output.append('mol : brush : freedom : restricted')
        output.append('mol : brush : theta : '+str(self.n*self.f*self.sigma))

        output.append('lat : mylat : geometry : flat')
        output.append('lat : mylat : n_layers : '+str(Zbox))
        output.append('lat : mylat : lambda : 0.166666666666666666666')
        output.append('lat : mylat : lowerbound : surface')
        output.append('mon : S : freedom : frozen')
        output.append('mon : S : frozen_range : lowerbound')
        output.append('/////////////////////////////////////////////')

        output.append('output : filename.out : type : ana')
        output.append('output : filename.out : write_bounds : false')

        output.append('/////////////////////////////////////////////')
        output.append('sys : noname : overflow_protection : true')
        output.append('newton : isac : iterationlimit : 100000')
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

    def plotProfile(self, profile_name = 'sys : noname : potential'):
        global g, xy, y_axis
        markerSize = '1pt'
        dataprefix = self.fname+'_'
        #~ dataprefix = ''
        try:
            type (g) == veusz.Embedded
            page = g.Root.page1
            graph0 = g.Root.page1.grid1.graph1
            graph1 = g.Root.page1.grid1.graph2
            graph2 = g.Root.page1.grid1.graph3
        except (NameError, AttributeError):
            g = veusz.Embedded(self.fname)
            g.EnableToolbar()

            page = g.Root.Add('page')

            grid = page.Add('grid')
            grid.bottomMargin.val ='0cm'
            grid.leftMargin.val = '0cm'
            grid.rightMargin.val = '0cm'
            grid.topMargin.val = '0cm'
            grid.rows.val = 2
            grid.columns.val = 1
            graph0 = grid.Add('graph')
            graph1 = grid.Add('graph')
            graph2 = grid.Add('graph')


        #~ g.To('/page1/grid1/graph1')

        self.loadData()
        #~ r=arange(1,len(sites))
        z = arange(0,self.Zbox)
        z = z[1:]
        x_dataname = dataprefix+'z'
        g.SetData(x_dataname, z)

        graph = graph0
        ########################
        # Plot polymer dencity #
        what_toPlot = 'phi'
        profile_toPlot = array(self.datadict['mol : brush : phi'])
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = 'black'
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'



        ####################
        # Plot counterions #

        what_toPlot = 'phi_Na'
        profile_toPlot = array(self.datadict['mol : na : phi'])
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = 'black'
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = 'red'
        ###############
        # Plot coions #

        what_toPlot = 'phi_Cl'
        profile_toPlot = array(self.datadict['mol : cl : phi'])
        #profile_toPlot=profile_toPlot[1:-1]
        y_dataname = dataprefix+what_toPlot

        g.SetData(y_dataname, profile_toPlot)
        xy = graph.Add('xy')
        xy.xData.val = x_dataname
        xy.yData.val = y_dataname
        xy.marker.val = 'none'

        xy.MarkerFill.color.val = 'black'
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'
        xy.PlotLine.color.val = 'blue'

        y_axis=graph.y
        x_axis=graph.x

        #y_axis.log.val = True
        #y_axis.min.val = 1e-4
        x_axis.max.val = 200
        y_axis.label.val = '\\phi_{Na}, \\phi_{Cl}, \\phi'
        x_axis.label.val = '\italic{z}'


        # Plot the distribution of the ends
        graph = graph1
        ########################
        # Plot polymer dencity #
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

        xy.MarkerFill.color.val = 'black'
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'


        y_axis=graph.y
        x_axis=graph.x

        #y_axis.log.val = True
        #y_axis.min.val = 1e-8
        x_axis.max.val = 200
        y_axis.label.val = '\\phi_{ends}'
        x_axis.label.val = '\italic{z}'




        # Plot the distribution of the branching points
        graph = graph2
        ########################
        # Plot polymer dencity #
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

        xy.MarkerFill.color.val = 'black'
        xy.markerSize.val = markerSize
        xy.ErrorBarLine.width.val = '2pt'
        xy.ErrorBarLine.transparency.val = 50
        xy.PlotLine.width.val = '1.5pt'
        xy.PlotLine.style.val = 'solid'


        y_axis=graph.y
        x_axis=graph.x

        #y_axis.log.val = True
        #y_axis.min.val = 1e-8
        x_axis.max.val = 200
        y_axis.label.val = '\\phi_{bp}'
        x_axis.label.val = '\italic{z}'



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
