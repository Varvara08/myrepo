import os
import numpy as np
#execfile('brush_l_class.py'); for line brush
execfile('brush_potention_class.py')
<<<<<<< HEAD
#mypath = '/home/alexander/imc/StarBrush/ph-sensitive/flatbrush/'
mypath = '/Users/varvara/Documents/IMC/ph-sensitive/' 
#mylist = [220] #np.arange(200,550,50)
mylist = np.arange(0.001,0.3,0.01)
=======
mypath = '/home/alexander/imc/StarBrush/ph-sensitive/flatbrush/'
mylist = np.arange(150,500,50)
>>>>>>> origin/master
#mylist = np.arange(0.25,0.8,0.25)
#mylist2= np.arange(1,2,1)

#####Dict and lambda
varpardictbp = {}; varpardictne = {}; varparfirstbp = {}; varparfirstne = {}; num_lambda = 0.166666666666666667;line_color_num = 1; # 5
line_type_num = 1; mylistav=[]; mylistdis=[];

##### Chain calc/sfbox options in potention###

chain_length = 220 

##### Star calc options in potention###

f1 = 3
arm_length = 100

##### Caclulation formula potention ###

n_formula = 100
f_formula = 3 # Its uses M formula
h2n = 0.25

##### Brush options ###

<<<<<<< HEAD
f = 3 # phi ---> potention ---> mod_pot ---> gt and gf.
n = 100
pK = 0.0 # it's was 5.0
pH = 7.0
val = 0.0 # it's was -1.0
cl = 1e-4
#sigma = 0.01 #
=======
f = 2 # phi ---> potention ---> mod_pot ---> gt and gf.
n = 100
pK = 0.0 # it's was 5.0
pH = 7.0
val = -0.3 # it's was -1.0
cl = 1e-4
sigma = 0.1 #
>>>>>>> origin/master

#####################
for vary_parameter in mylist:
#    for vary_parameter2 in mylist2:
 print '######### Brush option ####'
 print 'Forming a brush... \n f = '+str(f)+'\n n = '+str(n) + '\n pK = '+str(pK)+'\n val = '+str(val)+'\n sigma = '+str(vary_parameter)+'\n pH = ' + str(pH)
 print '######### Chain option #### \n chain_length = '+str(chain_length)
 print '######### Calc option #### \n f1 = '+str(f1)+'\n arm_length = ' + str(arm_length)
 print '######### Potention option #### \n h/2n = '+str(h2n)
 print '#########'
 b = brush()
 b.f = f; b.pK = pK ; b.val = val; b.n = n ; b.cl = cl ; b.sigma = vary_parameter ; b.pH= pH ; b.chain_length = chain_length ; b.h2n = h2n
 b.iguess = ''
 print 'SFBox time...'
 b.sfbox()
 b.chainInPotention()
 b.brushInPotention()
 #b.mylambda()
 b.calc_av_dis()
 print 'Plotting...'
 #b.plotProfile()
 #b.plotDiagram()
 b.plotAvDis()
 #g.Save(mypath+'data/'+b.fname+'.vsz')
# os.chdir('/home/alexander/mygit/UP/python')
 os.chdir('/Users/varvara/mygit/UP.eps/python')
 line_color_num+=1
