import os
import numpy as np
#execfile('brush_l_class.py'); for line brush
execfile('brush_d_class.py')
#mypath = '/home/alexander/imc/StarBrush/ph-sensitive/flatbrush/'
mypath = '/Users/varvara/Documents/IMC/pH-sensitive/'
mylist = [1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01]
#mylist2= np.arange(0,3.6,0.4)
#####Dict and lambda
varpardictbp = {}; varpardictne = {}; varparfirstbp = {}; varparfirstne = {}; line_color_num = 1; num_lambda =0.16666666666666666667;
#####

#chain_length = 100

#####Brush options###
f = 3
pK = 5.0
pH = 7.0
val = -1.0
n = 100
#cl = varya
sigma = 0.01
chi = 0.0
#####################
for vary_parameter in mylist:
    #for vary_parameter2 in mylist2:
     cant_solve_problem = 'False'
     print 'Forming a brush... \n f = '+str(f)+'\n n = '+str(n)+'\n pK = '+str(pK)+'\n val = '+str(val)+'\n sigma = '+str(sigma)+'\n pH = ' + str(pH) + '\n cl=' +str(vary_parameter) + '\n chi=' +str(chi)
     b = brush()
     b.f = f; b.pK = pK ; b.val = val; b.n = n ; b.sigma = sigma; b.pH= pH; b.cl=vary_parameter; b.chi=chi
     b.iguess = ''
     print 'SFBox time...'
#     if (b.pH == 9.5) and (b.sigma == 0.016) or (b.pH == 10.25) and (b.sigma == 0.016) or (b.pH == 13.0) and (b.sigma == 0.016) or (b.pH == 8.5) and (b.sigma == 0.019):
#	print 'There is bad point'
#     else :
     b.sfbox()
     if cant_solve_problem == 'True':
	print '#####'
	print 'It`s seems that SFBox does not find a solution'
	print '#####'
	print '\n WARNING: ALWAYSE CHECK FROM WHERE COMES DATA!!! '
     else:
        #b.chainInPotention()
        b.mylambda()
        print 'Plotting...'
        b.plotProfile()
        #b.plotDiagram()
        g.Save(mypath+'data/'+b.fname+'.vsz')
        os.chdir('/Users/varvara/mygit/UP.eps/python')
        line_color_num+=1
