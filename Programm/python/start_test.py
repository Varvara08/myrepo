import os
import numpy as np
#execfile('brush_l_class.py'); for line brush
execfile('brush_d_class.py')
mypath = '/Users/varvara/Documents/IMC/pH_brush/'
mylist=[1.0]
for a in mylist:
 b = brush()
 b.pK = a
 b.chi = 0.5
 b.iguess = ''
 b.val = -1
 b.sfbox()
 b.plotProfile()
 g.Save(mypath+'data/'+b.fname+'.vsz')
 	#~ y_axis.min.val = 1e-7
 	#~ del g
# os.chdir('/home/alexander/mygit/UP/python')