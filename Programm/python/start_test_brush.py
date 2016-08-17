import os
import numpy as np
#execfile('brush_l_class.py'); for line brush
execfile('brush_d+l_class.py')
mypath = '/home/alexander/imc/StarBrush/ph-sensitive/flatbrush/'
mylist = np.arange(0.1,0.4,10)
#mylist2= np.arange(1,2,1)

#####Dict and lambda
varpardictbp = {}; varpardictne = {}; varparfirstbp = {}; varparfirstne = {}; line_color_num = 1; num_lambda = 0.166666666666666667;

##### Chain options ###

chain_length = 250

##### Calc options ###

f1 = 3
arm_length = 100

##### Caclulation formula potention ###

n_formula = 100
f_formula = 3
h2n = 0.75

##### Brush options ###

f = 3
n = 100
pK = 0.0 # it's was 5.0
pH = 0.0
val = 0.0 # it's was -1.0
cl = 1e-1
#sigma = 0.095 #VARYA

#####################
for vary_parameter in mylist:
#    for vary_parameter2 in mylist2:
    print '######### Brush option ####'
    print 'Forming a brush... \n f = '+str(f)+'\n n = '+str(n) + '\n pK = '+str(pK)+'\n val = '+str(val)+'\n sigma = '+str(vary_parameter)+'\n pH = ' + str(pH)
    print '######### Chain option #### \n chain_length = '+str(chain_length)
    print '######### Calc option #### \n f1 = '+str(f1)+'\n arm_length = ' + str(arm_length)
    print '#########'
    b = brush()
    b.f = f; b.pK = pK ; b.val = val; b.n = n ; b.cl = cl ; b.sigma = vary_parameter ; b.pH= pH ; b.chain_length = chain_length
    b.iguess = ''
    print 'SFBox time...'
    b.sfbox()
    b.chainInPotention()
    b.brushInPotention()
    b.mylambda()
    print 'Plotting...'
    b.plotProfile()
    #b.plotDiagram()
    #g.Save(mypath+'data/'+b.fname+'.vsz')
    os.chdir('/home/alexander/mygit/UP/python')
    line_color_num+=1
