import os
import numpy as np
#execfile('brush_l_class.py'); for line brush
execfile('brush_d_class.py')
mypath = '/home/alexander/imc/StarBrush/ph-sensitive/flatbrush/'
mylist = np.arange(2,12,2)
#mylist2= np.linspace(4,9,20)
#####Dict and lambda###
varpardictbp = {}; varpardictne = {}; varparfirstbp = {}; varparfirstne = {}; line_color_num = 1; num_lambda =0.16666666666666666667;
#####ChainInPotention options###

chain_length = 100

#####Brush options###
f = 3
#pK = varya
pH = 7.0
val = -1.0
n = 100
cl = 1e-1
sigma = 0.01
#####################
for vary_parameter in mylist:
    #for vary_parameter2 in mylist2:
     print 'Forming a brush... \n f = '+str(f)+'\n n = '+str(n)+'\n pK = '+str(vary_parameter)+'\n val = '+str(val)+'\n sigma = '+str(sigma)+'\n pH = ' + str(pH)
     b = brush()
     b.f = f; b.pK = vary_parameter ; b.val = val; b.n = n ; b.cl = cl ; b.sigma = sigma ; b.pH= pH
     b.iguess = ''
     print 'SFBox time...'
     b.sfbox()
#     b.chainInPotention()
     b.mylambda()
     print 'Plotting...'
     b.plotProfile()
#     b.plotDiagram()
     g.Save(mypath+'data/'+b.fname+'.vsz')
     os.chdir('/home/alexander/mygit/UP/python')
     line_color_num+=1
