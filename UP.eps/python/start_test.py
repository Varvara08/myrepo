import os
import numpy as np
#execfile('brush_l_class.py'); for line brush
execfile('brush_d_class.py')
mypath = '/home/alexander/imc/StarBrush/ph-sensitive/flatbrush/'
mylist=[1.0,2.0,5.0,7.0,8.0,9.0]
for a in np.arange(0.,11.0):
 b = brush()
 b.pK = a
 b.chi = 0.5
 b.val = -1
 b.sfbox()
 b.plotProfile()
 g.Save(mypath+'data/'+b.fname+'.vsz')
 	#~ y_axis.min.val = 1e-7
 	#~ del g
# os.chdir('/home/alexander/mygit/UP/python')
