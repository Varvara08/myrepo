import os
import numpy
#execfile('brush_l_class.py');
execfile('brush_d_class.py');
execfile('functions.py');
b = brush()
mypath = '/home/alexander/imc/StarBrush/ph-sensitive/flatbrush/'
b.n = 100
b.chi = 0.0
#b.pK = 0.0

#for cl in logspace(log10(1e-6),log10(1e-2),10):
b.cl = 1e-1
b.sigma = 0.1
#b.val = -0.1
for pK in arange(0.,10.):
 b.sfbox()
 b.plotProfile()
 g.Save(mypath+'data/'+b.fname+'.vsz')
 	#~ y_axis.min.val = 1e-7
 	#~ del g
 os.chdir('/home/alexander/mygit/UP/python')
