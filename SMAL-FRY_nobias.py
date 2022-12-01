#############
# SMAL-FRY  #  __/|/|
############# (._>/\|  
#
# - SiMple AppLication of FisheR machinerY
#
# Developed by William Matthewson (Unige - william.matthewson@unige.ch) and Dennis Stock (Unige), based on a prior code from Mario Ballardini (University of Bologna).
# Used for simple Fisher forecasting and calculation of redshift-weighted number counts spectra in 2203.07414.
#
# The Parameters in this file are for a minimal example that can be used to quickly test that the various pieces of the code are working properly.
# Please see the following brief introduction for more details of this minimal example and how to use this code.
#
#############
# HOW TO USE:
#############
#
# THIS IS THE PARALLELISED VERSION AND REQUIRED ADDITIONAL EXPLANATION!!!
#
# The action of this code can be divided into three main sections:
# (1) - A modified version of CLASSgal is used to generate the spectra required. I.e. Galaxy number counts spectra (also referred to as \Delta). 
# (2) - A subroutine in this python code uses the generated \Delta spectra to produce the redshift-weighted number count spectra (also referred to as \zeta), which are then saved together with the \Delta spectra in new output files.  
# (3) - The full complement of spectra are read in and combined with user-defined parameters for observational surveys to produce Fisher forecasts for a set of parameters, also user-defined. Note: If adding new parameters for the Fisher forecasting, it is important that these parameters be linked somehow to the CLASS input, so that they effect the generated spectra and subsequent numerical derivatives used in the forecasts.
# A new user will need to familiarise themselves with the input parameters of CLASS, the survey (biases, galaxy distributions, sky coverage etc.) they would like to forecast and the available parameters that are supported by this code. For now the code has only been tested with spectroscopic galaxy surveys, but the option exists to generalize to photometric/Intensity mapping surveys in the Fisher code and the modified version of CLASSgal. 
# The major modifications to CLASSgal are: 
# - The addition of certain window functions (see \omega(z_i,z) in the paper 2203.07414) that may be used to generate the \zeta spectra. 
# - The addition of a 'selection_list' variable that allows the window function for each bin to be specified individually.
# - The inclusion of redshift-dependent function forms for the linear bias, b(z), specific to certain surveys, eventually with the aim of marginalising over the fitting parameters for the biases. 
#
# One important step, which is how the user controls the chosen survey, is to set the elements of 'N' to 0, EXCEPT for the survey of interest, where it should be 2x the number of \Delta bins. This particular architectural choice is to allow for multitracing to be supported more easily at a later time.
#
# There is also a "template" file included, namely 'Deltazeta_template.ini', where any remaining CLASS parameters are set that do not need to change between the different surveys.
#
# Once the survey specifications have been chosen, the user can run (individually or in one go*) via this python code: the modified CLASS code and the \zeta spectra extraction. 
# *- By choosing the values of 'RUN_eSPECTRA' and 'EXTRACT_ZETA', respectively. It may make sense at the beginning to run the CLASS part alone at first for a smaller run to make sure that accuracy parameters etc. are ideal and the spectra being generated are what they are required to be.
# The CLASS part will generate \Delta spectra for each bin, using 3 different window functions: 2 for the extraction of the \zeta spectra, and one as the \Delta spectrum with a Gaussian window function. 
# This means that the eventual 'extracted' spectra will be \Delta_i\Delta_j, \zeta_i\zeta_j, \Delta_i\zeta_j and \zeta_i\Delta_j for (i=j) auto- and (i!=j) cross-correlations of the chosen set of bins centred at redshifts z_i. (I.e. the SAME bins will be used for \Delta and \zeta spectra, which is also the reason why N contains 2x the number of bins for an individual spectra.) 
#
# Finally, if all the required spectra have been generated, the code should complete by performing the Fisher forecast. It is recommended that the user check the numerical derivatives ('DERIVATIVE CHECK') to make sure that the chosen step-size ('step') is stable numerically.
#
# In its current state, the code will output the Fisher forecast constraints and covariance matrices for:
# - The full \Delta\&\zeta information
# - The \Delta spectra alone
# - The \zeta spectra alone
# - The combination of \Delta\Delta and \zeta\zeta correlations (i.e. the first output on this list, without the cross \Delta\zeta correlations).
#
# A simple summary for a run:
# - Download the modified CLASSgal code directory, this python file.
# - Change 'FILEPATHTO' to the correct one that corresponds to the CLASS directory, don't forget to end the string with a slash '/'.
# - [make clean], [make class] in the CLASS directory so that the modified version is compiled.
# - Within this directory, if necessary, create the structure for outputs, see 'FILE STRUCTURES AND LABELLING'.
# - Run the minimal example to check everything is in place.
# - [This may be included in the next steps]: Undo the 'damage' done to set the code up for the minimimal example, search 'undo the minimal example'.
# - Choose a survey and enter all the relevant parameters into this python file.
# - Check to see that the Fisher parameters are correct.
# - run [python3 SMAL-FRY.py] with 'RUN_SPECTRA' and 'EXTRACT_ZETA' both set equal to 1.
# - Depending on the spectra concerned, accuracy, number of bins, type of numerical derivatives etc. you may have some time on your hands at this stage.
# - Revel in your newly-completed constraints.
# - Once finished the previous step, make sure to check that the CLASS spectra look realistic, the derivatives are stable and that all the spectra are (a) the ones you think they are (check 'NO') and (b) are accurate and correct. 
# - If you complete the previous step without issue, you may return to enjoying your constraints.
#
################################
# FILE STRUCTURES AND LABELLING:
################################
#
# The user should download the modified version of CLASSgal from the same GitHub page (https://github.com/WillMatt4/SMAL-FRY), and compile it. Search 'FILEPATHTO' in this code to replace it with the relevant directory address. 
# If it does not exist, you will need to create inside the directory 'class_public-3.0.1_mod' a directory structure like: '/spectra/'+suptype, where suptype is replaced by the string value of 'suptype'.
# Each time you run CLASS, it should produce in this sub-directory the \Delta cls data files, suffixed with a two digit number, starting at 00, and prefixed with the label 'suptype'. It is important that you ensure that the value of 'NO' correctly corresponds to the CLASS \Delta spectra that you are basing the analysis on.
#
#
##################
# MINIMAL EXAMPLE:
##################
#
# In the current state, the code will produce a 2 bin forecast for ONLY 3 cosmological parameters, from Euclid, with the bins evenly dividing the redshift range of [0.9, 1.8]. The maximum ell computed is also set only to 300, with the redshift-dependent lmax cut-off commented out. This is by no means intended to produce good constraints, but should provide evidence of whether the code has been correctly interpreted and set-up, or where there might be problems.
#
# In order to undo the minimal example and set up for your first real run, you will need to:
# - Choose the number of bins in your survey and set in the relevant element of 'N'.
# - Make sure that all the parameters you want are included in the parameter dictionary, 'params = {}'
# - Make sure that these parameters are not overwritten, and instead get dynamically input into CLASS. Search 'Minimal example background parameters' and comment out the 2 lines relevant to the minimal example, and uncomment the 2 lines relevant to 'Otherwise'.
# - Set 'lmax' to a value above the largest z-dependent one required for your bins (or even a bit larger, for CLASS accuracy reasons.)
# - Search 'z-dependent lmax' in this python file and turn this feature on as described there.
# - Delete the generated CLASS output files (suffix '00') or remember to change 'NO' to '01' for the next run.
#
# Certain features that are envisioned, but not yet included are tagged with '[TBC]' and a simpmle search of this python file will reveal these to the user.
#
from numpy import *
import numpy as np
#from pylab import *
from scipy.special import erf
from scipy.integrate import quad
from scipy import interpolate
import time
import os.path
from os import path
from os import system
#import sh as pbs
import mpi4pywrapper as mpi
start=time.time()
#mpi.barrier = lambda : 1


##################################################################
##                                                              ## 
## N.B.!  WHEN CHOOSING SURVEY, CHANGE VALUES OF 'NO' AND 'N'   ##
##                                                              ##
##################################################################

################################
#COSMOLOGY AND INPUTS FOR CLASS:
################################
RUN_SPECTRA = 0                                   # Should CLASS be used to generate the necessary spectra? YES: 1, NO: 0.
EXTRACT_ZETA = 0                                  # Should the redshift-weighted number count power spectrum? YES: 1, NO: 0.
RUN_FISHER = 1                                    # Should the Fisher analysis be run (DOESN'T REQUIRE MULTIPLE CORES)
CLASSPATH = '/home/users/m/matthews/scratch/Zeta/SMAL-FRY/class_public-3.0.1_mod/' #
OUTPATH = '/home/users/m/matthews/scratch/Zeta/outputs/' #
NO = '06'                                         # This Number relates to the filename and numbering scheme of the input CLASS spectra.
 
#################
#INITIALIZATIONS:
#################

#Cosmological background parameters:
####################################
params={}
params['wb'] = 0.022383                          # \Omega_b*h^2
params['wc'] = 0.12011                           # \Omega_c*h^2
params['h'] = 0.6732
params['ns'] = 0.96605
params['logAs'] = 3.0448                        # log(10^{10}A_s
k_pivot = 0.05                                   # [h/Mpc]
tau0 = 14187.020887                              # Age of universe [Mpc].
lmax = 1500
kmax = 0.2                                       # Wave number [Mpc^-1] relating to the non-linear cut-off in ell for each bin used in the Fisher analysis.


#Fisher and Survey parameters:  
##############################

# The way that the code deals with different surveys is that the user inputs all the relevant parameters which can be stored and then accessed by changing the values of the variables "NO" and "N". The former controls the CLASS spectra that are read-in and the latter controls the chosen survey. More specifically, the user sets the parameters for the nth survey in the nth element of the various arrays in this section, controlling which survey is used in the forecast by setting the corresponding element for the number of bins in "N" to a value>0, and the rest to 0.
# Currently, with the addition of the redshift-weighted number count power spectrum, the code does not completely support multi-tracer forecasts [TBC], though this generalisation should be fairly straight-forward. 
#

suptype = 'Deltazeta'                            # Label for CLASS and Fisher output files.
ctype = ['gal_spec','gal_spec']                  # Flag that controls the type of survey, spectroscopic galaxies, photometric galaxies, HI intensity mapping for the purposes of CLASS spectra generation. 
subctype = ['Euc_spec','SKA_2_spec']             # Flag that controls the specific survey biases, dNdz etc. to be used. Used as a label in the output files of Fisher forecasts. 
zMx = [1.8,2.0]                                  # Maximum redshift bound for each survey
zMn = [0.9,0.1]                                  # Minimum redshift bound for each survey
N = [40,0]                                        # 2x the individual number of bins for the survey considered. All other elements should be set to zero, unless the code has been adapted for multi-tracer calculations [TBC]. Since each bin will have a number counts power spectrum (\Delta) and a redshift-weighted number counts power spectrum (\zeta), this should always be an EVEN number = 2x Number of C_\ell^{\Delta} bins!
Nsum = sum(N)
Redsig = [0.001,0.001]                           # Redshift accuracy, later should be multiplied (1+z) for each bin. Redshift window widths, on the other hand, are stored in zsigmas [TBC] #CURRENTLY this parameter has no part in the calculation, but e.g. for photometric galaxies where the resolution is more limiting, should be taken into account in the galaxy distribution function.
#sigmaZ = Redsig[0]
biastype = [2,3]                                 # Linear bias general function form: 0 - Constant bias: b=biasno 1 - DESI bias: b(z) = biasno/D(z) 2 - Euclid bias: b(z) = biasno + biasno2*z 3 - SKAII bias: b(z) = biasno*exp(biasno2*z).
biasno = [0.79,0.5887] 
biasno2 = [0.68,0.8130]
surv = [i for i in range(len(N)) if N[i]>0]

#params['b1'] = biasno[surv[0]]                         #1st bias parameter
#params['b2'] = biasno2[surv[0]]                        #2nd bias parameter
#print('\n\nb1: '+str(params['b1'])+' b2: '+str(params['b2'])+'\n\n')


#Noise and Fisher limits:
#########################

# NB! ONLY MAXIMUM OF ONE (Marg OR Fix) MAY BE NON-ZERO:
# The Fix and Marg parameters have NOT been tested in this version of the code and should be checked if used [TBC]. 
Fix = 0                                          # Number of parameters to fix (remove before inversion) (starting at the parameter with the highest index).
Marg = 2                                         # Number of parameters to marginalise (remove after inversion) over (starting at the parameter with the highest index)
# Currently, if left as default (0), the code will produce the conditional and marginalized errors (but no additional parameters are marginalized over).
Cond = 1                                         # Should conditional errors be found? YES: 1, NO: 0.


Noise = 1                                        # Should shot(/instrumental in some surveys) Noise be included? YES: 1, NO: 0.
Lmin_TrPre = [10,10]                             # A minimum ell bound for the Fisher analysis to account for foreground cleaning/other effects that mean the lowest ells cannot be used.
Lmin_Tr = []
for i in range(len(Lmin_TrPre)):
    Lmin_Tr += [Lmin_TrPre[i] for j in range(N[i])]

fskydegsq = [15000,30000]                        # The sky coverage of each survey in deg^2 
nParsLess = 5                                    #Number of parameters excluding parameters to be marginalised over. This part of the code is not yet active [TBC].

SpectraPath = OUTPATH+'spectra/'+suptype+'/'   # Path to the CLASS-generated spectra. If the earlier instructions about labelling are followed, this can be left unchanged. 
f=open(suptype+'_fiducial_params.dat','w')
for par in params:
   f.write(str(params[par])+'\n')
f.close()

#Redshift distributions (tracer/survey-dependent):
##################################################

#
# Here the various parameters relating to the redshift distributions can be set.
#

def dNdz(z,Tracer):
    if (ctype[Tracer]=='gal_phot' or ctype[Tracer]=='gal_spec'):
        if (subctype[Tracer]=='Euc_spec'):
            dndz = (-1494.0710866859654*z**2+2598.212251011541*z+717.0942492893472) #quadratic fit from 1910.09273 for galaxy number density [in deg^-2]

    elif (ctype[Tracer]=='HI'):
        dndz = 0.055919 + 0.23242*z - 0.024136*z**2
    return dndz


def get_nth_key(dictionary, n=0):
    if n < 0:
        n += len(dictionary)
    for i, key in enumerate(dictionary.keys()):
        if i == n:
            return key
    raise IndexError("dictionary index out of range") 

Keys = [get_nth_key(params,i) for i in range(len(params))]


#Numerical derivatives:
#######################
stencil= 5 #3 or 5 point stencil
step=0.01

WHICHTRACERS=[] #For use especially in later multi-tracer forecasts: a useful list containing the array indices for the survey(s) considered. 
for i in range(len(N)):
  if N[i] != 0: 
    WHICHTRACERS.append(i)

Tracer = WHICHTRACERS[0]
Lmax_Tr = [lmax for i in range(N[Tracer])]



#Derivative Functions:
######################
#p is the parameter name: ('wb'), ('wc'), ('H0'), ('ns'), ('logAs') etc.
#X is the spectra index: 1:(d1d1 ...)... Nbins:(d1dN...) Nbins+1:(d1z1...)... Nbins+N=2N:(d1zN...)  etc. 
def deriv(X,p,l):
    if stencil==3:
        return (Cl_p[p][l-2,X]-Cl_m[p][l-2,X])/(2*step*params[Keys[p]])
    elif stencil==5:
        return (-Cl_p2[p][l-2,X] + 8*Cl_p[p][l-2,X]-8*Cl_m[p][l-2,X]+Cl_m2[p][l-2,X])/(12*step*params[Keys[p]])
    
def derivTr(X,p,l):
    if stencil==3:
        return (Tr_p[p][l-2,X]-Tr_m[p][l-2,X])/(2*step*params[Keys[p]])
    elif stencil==5:
        return (-Tr_p2[p][l-2,X] + 8*Tr_p[p][l-2,X]-8*Tr_m[p][l-2,X]+Tr_m2[p][l-2,X])/(12*step*params[Keys[p]])


########
#HEADER:
########
print('\n---|----------------------------------------------------------------------------|-')
print('-|----------------------------------------------------------------------------|---')
print('---|  SSSSS   MMM    MM    AA    LL          FFFFFFFF RRRRRR   YY     Y       --|-')
print('-|-- S    SS  MM MM M M   AA A   LL          FF       RR   RR   YY   Y        |--')
print('---|  SS      MM  MM  M  AA   A  LL     ___  FF       RR    RR   YY Y         --|-')
print('-|--   SSS    MM      M AAAAAAAA LL    |___| FFFFF    RR   R      YY          |---')
print('---|     SS   MM      M AA     A LL          FF       RRRRR       YY          --|-')
print('-|-- SS    S  MM      M AA     A LL          FF       RR   RR     YY   __/|/| |---')
print('---|  SSSSS   MM      M AA     A LLLLLLLL    FF       RR     RR   YY  (._>/\| --|-')
print('-|----------------------------------------------------------------------------|---')
print('---|----------------------------------------------------------------------------|-')
print(' - SiMple AppLication of FisheR machinerY.')
print('\n\n')




#############
#CLASS INPUT:
#############

#
# These SKA(II) parameters are set manually, otherwise the survey binning will be performed across the redshift range in equal redshift portions, according to the number of bins set. The default redshift (Gaussian) bin width is set to the full-width at half maximum. NOTE: only gaussian window functions are fully-supported by the \zeta and Fisher code, even though CLASS handles other window types. [TBC]  
#
SKAmeans = [0.2, 0.28, 0.36, 0.44, 0.53, 0.64, 0.79, 1.03, 1.31, 1.58, 1.86]
SKAsigmas = [0.1, 0.08, 0.08, 0.08, 0.1, 0.12, 0.2, 0.28, 0.28, 0.28, 0.28]

#Window and Redshift function (including for \zeta):
###################################################
for Tracer in WHICHTRACERS:
      bins=[i for i in zeros(int(N[Tracer]/2))] 
      zmeans=[i for i in zeros(int(N[Tracer]/2))] 
      zsigmas=[i for i in zeros(int(N[Tracer]/2))] 
      for i in range(int(N[Tracer]/2)):
            bins[i]=[zMn[Tracer]+i*(zMx[Tracer]-zMn[Tracer])/(N[Tracer]/2),zMn[Tracer]+(i+1)*(zMx[Tracer]-zMn[Tracer])/(N[Tracer]/2)]
            zmeans[i] = mean(bins[i]) 
            zsigmas[i] = (bins[i][1] - bins[i][0])/2/np.sqrt(2*np.log(2)) 
                                                 #The bin half-width is half the full width at half maximum, thus for the gaussian window: sigma = (bin width)/2/sqrt(2*ln(2))

            #Manual bin parameters
            if(subctype[Tracer] == 'SKA_2_spec'):
                bins[i] = [SKAmeans[i]-SKAsigmas[i]/2,SKAmeans[i]+SKAsigmas[i]/2]
                zsigmas[i] = SKAsigmas[i]/2/np.sqrt(2*np.log(2))
                zmeans[i] = SKAmeans[i]

#Noise calculation initialisation:
##################################

NprevC = np.zeros(len(N))                        
for Tracer in WHICHTRACERS:
  NprevC[Tracer] = int(sum(N[0:Tracer]))
  fskyTr= fskydegsq[Tracer]/(4*pi*(180/pi)**2)

#########
#SPECTRA:
#########

#(We have to divide Nsum by 2 because we only want Nsum/2 of each bin: \Delta AND \zeta)
Nsumo2 = int(Nsum/2)
corrDic = []
finalheader = '# ell'
for i in range(Nsumo2):
  for t in range(2):
    for j in range(Nsumo2):
       if t==0:      #di x dj
          corrDic.append(' D_'+str(i+1)+'xD_'+str(j+1)) 
          finalheader += ' D_'+str(i+1)+'xD_'+str(j+1)
       elif t==1:    #di x zj
          corrDic.append(' D_'+str(i+1)+'xZ_'+str(j+1)) 
          finalheader += ' D_'+str(i+1)+'xZ_'+str(j+1)
for i in range(Nsumo2):
  for t in range(2):
    for j in range(Nsumo2):
       if t==0:      #zi x dj
          corrDic.append(' Z_'+str(i+1)+'xD_'+str(j+1)) 
          finalheader += ' Z_'+str(i+1)+'xD_'+str(j+1)
       elif t==1:    #zi x zj
          corrDic.append(' Z_'+str(i+1)+'xZ_'+str(j+1)) 
          finalheader += ' Z_'+str(i+1)+'xZ_'+str(j+1)


Tracer = WHICHTRACERS[0] # Fixed here because we only have one tracer, [TBC] for multi-tracer forecasts in the future.

def Spectra(steptype,params,key):
# This function attributes the kind of spectrum (fiducial, +/- 1*step, +/- 2*step to be generated, and calls RUNCLASS() to generate it.)
    params0=params.copy()
    if steptype<0:
       signV = 'minus'+str(abs(steptype))
    elif steptype==0:
       signV = 'fiducial'
    else:
       signV = 'plus'+str(abs(steptype))

    if steptype==0:
        RUNCLASS(signV,params0,key)
    else:
        params0[key]=params0[key]*(1+steptype*step)
        RUNCLASS(signV,params0,key) 
    #itrs = mpi.rank
    #if itrs>0:
      #print(itrs,params0[key],steptype,step)
def RUNCLASS(signV,params,key):
# This function is responsible for running the CLASS code to generate and save the necessary \Delta spectra.
    if signV=='fiducial':
       print('\n\nRunning CLASS spectra: fiducial...\n\n')
    else:
       print('\n\nRunning CLASS spectra: '+key+' '+signV+'...\n\n')


    selectionmeans = ''
    selectionbiastypes = ''
    selectionbiass = ''
    selectionbias2s = ''
    selectionlists = ''
    selectionwidths = ''
    selectionmagbiass = ''

    for i in range(len(bins)):
         
         # Linear bisa b(z) and magnification bias s(z) parameters for input to CLASS.
         BIASNO = biasno[surv[0]]# params['b1']#
         BIASNO2 = biasno2[surv[0]]#params['b2']#
         #print(BIASNO,BIASNO2)
         s0 = -0.106875
         s1 = 1.35999
         s2 = -0.620008
         s3 = 0.188594
         sSKAII = s0 + s1*zmeans[i] + s2*zmeans[i]**2 + s3*zmeans[i]**3
         sEUC = 2/5*(0.583 + 2.02*zmeans[i] -0.568*zmeans[i]**2 + 0.0411*zmeans[i]**3)
         if(subctype[Tracer] == 'SKA_2_spec'):
            ssurvey = sSKAII
         elif(subctype[Tracer] == 'Euc_spec'):
            ssurvey = sEUC

         # In order to generate the spectra necessary for each bin requires 2 \Delta spectra with different kinds of window functions, in addition to the Gaussian case.
         # Each redshift bin's spectrum has to be generated 3 times, with window functions 3 - \omega(z_i,z) (in the paper 2203.07414), 4 - z*\omega(z_i,z) (in the paper 2203.07414), 0 - Gaussian. The "selection_list" parameter passed to the CLASS .ini file has been added to the CLASS code, along with the additional window functions, to allow for the control of each bin's window function individually. This is the reason for the following rather tedious generation of strings.
      
         if i==len(bins)-1:
            selectionmeans+=str(zmeans[i])+', '+str(zmeans[i]+1e-7)+', '+str(zmeans[i]+2e-7)
            selectionbiastypes += str(biastype[Tracer])+', '+str(biastype[Tracer])+', '+str(biastype[Tracer])
            selectionbiass += str(BIASNO)+', '+str(BIASNO)+', '+str(BIASNO)
            selectionbias2s += str(BIASNO2)+', '+str(BIASNO2)+', '+str(BIASNO2)
            selectionlists += '3,4,0' 
            selectionwidths += str(zsigmas[i])+','+str(zsigmas[i])+','+str(zsigmas[i])
            selectionmagbiass += str(ssurvey)+','+str(ssurvey)+','+str(ssurvey)
         else: 
            selectionmeans+=str(zmeans[i])+', '+str(zmeans[i]+1e-7)+', '+str(zmeans[i]+2e-7)+', '
            selectionbiastypes += str(biastype[Tracer])+', '+str(biastype[Tracer])+', '+str(biastype[Tracer])+', '
            selectionbiass += str(BIASNO)+', '+str(BIASNO)+', '+str(BIASNO)+', '
            selectionbias2s += str(BIASNO2)+', '+str(BIASNO2)+', '+str(BIASNO2)+', '
            selectionlists += '3,4,0,'
            selectionwidths += str(zsigmas[i])+','+str(zsigmas[i])+','+str(zsigmas[i])+','
            selectionmagbiass += str(ssurvey)+','+str(ssurvey)+','+str(ssurvey)+','




    system('cp '+CLASSPATH+suptype+'_template.ini'+' '+CLASSPATH+suptype+'_'+key+'_'+signV+'.ini')
    f=open(CLASSPATH+suptype+'_'+key+'_'+signV+'.ini','a') 
    f.write('l_switch_limber_for_nc_local_over_z = '+str(300)+'\n')  #These accuracy parameters can be adjusted to attain the necessary accuracy at lower ells, especially for narrow windows where Limber is not effective on large scales.
    f.write('l_switch_limber_for_nc_los_over_z = '+str(300)+'\n')    #These accuracy parameters can be adjusted to attain the necessary accuracy at lower ells, especially for narrow windows where Limber is not effective on large scales.
    f.write('h = '+str(params['h'])+'\n') 
    f.write('omega_b = '+str(params['wb'])+'\n') 
    f.write('omega_cdm = '+str(params['wc'])+'\n') 

    #Minimal example background parameters:
    #f.write('ln10^{10}A_s = '+str(3.0448)+'\n') 
    #f.write('n_s = '+str(0.96605)+'\n') 
    #Otherwise:
    f.write('ln10^{10}A_s = '+str(params['logAs'])+'\n') 
    f.write('n_s = '+str(params['ns'])+'\n') 

    f.write('k_pivot = '+str(0.05)+'\n') 
    f.write('number count contributions = density, rsd, lensing, gr\n') 
    f.write('non_diagonal = '+str(3*N[Tracer]/2-1)+'\n') 
    f.write('selection_mean = '+selectionmeans+'\n')
    f.write('selection_bias_type = '+selectionbiastypes+'\n')
    f.write('selection_bias = '+selectionbiass+'\n')
    f.write('selection_bias2 = '+selectionbias2s+'\n')
    f.write('selection_list= '+selectionlists+'\n') 
    f.write('selection_width = '+selectionwidths+'\n') 
    f.write('l_max_lss = '+str(lmax)+'\n')   
    if(signV=='fiducial'):   
       f.write('write background = yes\n')      
       f.write('write parameters = yes\n')
    else:     
       f.write('write background = no\n')      
       f.write('write parameters = no\n')
    f.write('selection_magnification_bias = '+selectionmagbiass+'\n')
    f.write('root = '+SpectraPath+suptype+'_'+key+'_'+signV+'_'+NO+'\n')
    f.close()
    
    system('(cd '+CLASSPATH+' && exec ./class '+suptype+'_'+key+'_'+signV+'.ini)')
    print('\n'+key+' '+signV+' spectra complete.\n')


def ExtractClDz(key,signV,zbar):
# This function is used to extract \Delta, \zeta and cross spectra:

    data2 = np.loadtxt(SpectraPath+suptype+'_'+key+'_'+signV+'_'+NO+'00_cl.dat')
    N2do4 = int(N[Tracer]**2/4)
    Cl_fz_zi_fz_zj = [[0 for j in range(N2do4)] for i in range(N2do4)] #extract the element like C[i][j] 
    Cl_fz_zi_f_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]
    Cl_f_zi_fz_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]
    Cl_f_zi_f_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]

    Cl_f_zi_d_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]
    Cl_d_zi_f_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]
    Cl_fz_zi_d_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]
    Cl_d_zi_fz_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]

    Cl_di_dj = [[0 for j in range(N2do4)] for i in range(N2do4)]
    Cl_zi_dj = [[0 for j in range(N2do4)] for i in range(N2do4)]
    Cl_di_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]
    Cl_zi_zj = [[0 for j in range(N2do4)] for i in range(N2do4)]


    #The indices here are used to convert the CLASS format [bin1bin1,bin1bin2,...,bin1binN, bin2bin2,bin2bin3,...,bin2binN,... ..., binNbinN], where each bin (cross-)correlation appears only once, into "CAMB format", which goes through all indices, like a matrix: [bin1bin1,bin1bin2,...,bin1binN, bin2bin1,bin2bin2,bin2bin3,... ..., binNbinN].
    def CLASSindConvert(Ind,N):
      j = Ind-1-np.floor((Ind-1)/3/N)*3*N
      i = np.floor((Ind-1)/3/N)
      if i<=j:
          return int(3*N*i + i + 1 - i*(i+1)/2+(j - i))
      else:
          return  int(3*N*j + j + 1 - j*(j+1)/2+(i - j))
      
   
    def IndexFunc(i,j,A,B,N):
      return CLASSindConvert(3*N*(3*(i)+A) + 3*(j)+B,N)
    


    # Here we generate the necessary precursors for all the \Delta, \zeta and cross spectra:
    N0o2 = int(N[Tracer]/2)
    for i in range(N0o2):
       for j in range(N0o2):
         Cl_fz_zi_fz_zj[i][j] =  data2[:, IndexFunc(i,j,1,2,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi
         Cl_fz_zi_f_zj[i][j] = data2[:, IndexFunc(i,j,1,1,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi
         Cl_f_zi_fz_zj[i][j] = data2[:, IndexFunc(i,j,0,2,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi
         Cl_f_zi_f_zj[i][j] = data2[:, IndexFunc(i,j,0,1,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi

         Cl_fz_zi_d_zj[i][j] = data2[:, IndexFunc(i,j,1,3,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi
         Cl_f_zi_d_zj[i][j] = data2[:, IndexFunc(i,j,0,3,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi
         Cl_d_zi_fz_zj[i][j] = data2[:, IndexFunc(i,j,2,2,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi
         Cl_d_zi_f_zj[i][j] = data2[:, IndexFunc(i,j,2,1,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi

         Cl_di_dj[i][j] = data2[:, IndexFunc(i,j,2,3,N0o2)]/(data2[:, 0]*(data2[:, 0]+1))*2*np.pi

    for i in range(N0o2):
       for j in range(N0o2):
         Cl_zi_dj[i][j] = zbar[i]*(Cl_fz_zi_d_zj[i][j]-Cl_f_zi_d_zj[i][j])
         Cl_di_zj[i][j] = zbar[j]*(Cl_d_zi_fz_zj[i][j]-Cl_d_zi_f_zj[i][j])
         Cl_zi_zj[i][j] = zbar[i]*zbar[j]*(Cl_fz_zi_fz_zj[i][j]-Cl_fz_zi_f_zj[i][j]-Cl_f_zi_fz_zj[i][j]+Cl_f_zi_f_zj[i][j])
    return Cl_di_dj,Cl_di_zj,Cl_zi_dj,Cl_zi_zj

def SaveClDz(key,signV,Cl_di_dj,Cl_di_zj,Cl_zi_dj,Cl_zi_zj):
# This function saves the extracted spectra to the final data files to be used in the Fisher analysis.

# As mentioned above, the format of this is like CAMB (NOTE: there are N \Delta (di) bins and N \zeta (zi) bins to be taken account of!)
# [Header]
#   #Cl_d1_d1 Cl_d1_d2 ... Cl_d1_dN Cl_d1_z1 ... Cl_d1_zN Cl_d2_d1 Cl_d2_d2 ... Cl_d2_dN Cl_d2_z1 ... Cl_d2_zN Cl_dN_d1 Cl_dN_d2 ... Cl_dN_dN Cl_dN_z1 ... Cl_dN_zN Cl_z1_d1 Cl_z1_d2 ... Cl_z1_dN Cl_z1_z1 ... Cl_z1_zN Cl_zN_d1 Cl_zN_d2 ... Cl_zN_dN Cl_zN_z1 ... Cl_zN_zN
# [Data]
# l = 2
#   .
#   .
#   .
# l = lmax

    Nsumo2 = int(Nsum/2)
    f=open(SpectraPath+suptype+'_'+key+'_'+signV+'_'+NO+'00_cl_Delta_zeta.dat','w')
    f.write(finalheader+'\n')
    for ell in range(2,lmax+1):
       f.write(str(ell)+' ')
       for i in range(Nsumo2):
         for t in range(2):
             for j in range(Nsumo2):
                if t==0:      #di x dj
                   f.write(str(Cl_di_dj[i][j][ell-2])+' ')
                elif t==1:    #di x zj
                   f.write(str(Cl_di_zj[i][j][ell-2])+' ')
       for i in range(Nsumo2):
         for t in range(2):
             for j in range(Nsumo2):
                if t==0:      #zi x dj
                   f.write(str(Cl_zi_dj[i][j][ell-2])+' ')
                elif t==1:    #zi x zj
                   f.write(str(Cl_zi_zj[i][j][ell-2])+' ')


      
       f.write('\n')



#Calculate Initial Spectra:
############################
if RUN_SPECTRA==1:
    
  #FIDUCIAL:
  ##########
  itr = mpi.rank
  if itr==0:
      Spectra(0,params,'')
  
  #mpi.barrier #wait for fiducial run to complete (background to be generated)
    
  if itr!=0:  
  #if itr>=0:  
      key = Keys[itr-1]
  #    key = Keys[itr]
  #PLUS STEP:
  ###########
      Spectra(1,params,key)
  #MINUS STEP:
  ############
      Spectra(-1,params,key)

      if stencil==5:

    #PLUS 2*STEP:
    #############
        Spectra(2,params,key)
    #MINUS 2*STEP:
    ##############
        Spectra(-2,params,key)
      
  mpi.barrier

Z,Ch,Hub,D,f  = np.genfromtxt(SpectraPath+suptype+'__'+'fiducial_'+NO+'00_background.dat',unpack=True,usecols = [0,4,3,17,18])

Dofz = interpolate.interp1d(Z, D, kind='cubic',axis = 0)    #Mpc
fofz = interpolate.interp1d(Z, f, kind='cubic',axis = 0)    #Mpc
Chif = interpolate.interp1d(Z, Ch, kind='cubic',axis = 0)   #Mpc
Hubf = interpolate.interp1d(Z, Hub, kind='cubic',axis = 0)  #/Mpc

zs = [mean(bins[i]) for i in range(len(bins))]
zbars = np.zeros(len(zs))
for i in range(len(zs)):
     fzi = lambda z: np.exp(-1/2*(z-zs[i])**2/zsigmas[i]**2)/np.sqrt(2*np.pi)/zsigmas[i]*z*Chif(z)**2/Hubf(z)
     fi=  lambda z: np.exp(-1/2*(z-zs[i])**2/zsigmas[i]**2)/np.sqrt(2*np.pi)/zsigmas[i]*Chif(z)**2/Hubf(z)
     # The factor of 30*zsigmas is arbitrary, it must be selected based on the characteristics of the survey and its window functions, so that the numerical integral is not over to large a range that may result in the resolution of the non-zero values of the amplitude being poor. This is the integral note.
     zbars[i]= quad(fzi, max([0,zs[i]-30*zsigmas[i]]), zs[i]+30*zsigmas[i])[0]/quad(fi, max([0,zs[i]-30*zsigmas[i]]), zs[i]+30*zsigmas[i])[0]



#Read and Initialize background
################################
if EXTRACT_ZETA==1:  
  
   ##key = Keys[7-1]
   ##for signV in ['minus1','plus1','minus2','plus2']:
   ##   Cl_di_dj,Cl_di_zj,Cl_zi_dj,Cl_zi_zj=ExtractClDz(key,signV,zbars)
   ##   SaveClDz(key,signV,Cl_di_dj,Cl_di_zj,Cl_zi_dj,Cl_zi_zj)

   

  itr = mpi.rank
  if itr==0:  
    Cl_di_dj,Cl_di_zj,Cl_zi_dj,Cl_zi_zj = ExtractClDz('','fiducial',zbars)
    SaveClDz('','fiducial',Cl_di_dj,Cl_di_zj,Cl_zi_dj,Cl_zi_zj)

  if itr!=0:   
   key = Keys[itr-1]
   for signV in ['minus1','plus1','minus2','plus2']:
      Cl_di_dj,Cl_di_zj,Cl_zi_dj,Cl_zi_zj=ExtractClDz(key,signV,zbars)
      SaveClDz(key,signV,Cl_di_dj,Cl_di_zj,Cl_zi_dj,Cl_zi_zj)

   mpi.barrier

###############################
#LET THE FISHERING COMMENCE!!!#
###############################

if RUN_FISHER==1:

#Check for missing spectra
############################ 
# This check might give false hope if the user forgets to set 'NO' to the correct one. 
  filedir = SpectraPath
  print('Reading from: '+filedir)

  nono=0
  if(path.exists(filedir+suptype+'__fiducial_'+NO+'00_cl_Delta_zeta.dat')==False):
    print('Missing fiducial spectra\n')

  for par in params:
    if(path.exists(filedir+suptype+'_'+par+'_minus1_'+NO+'00_cl_Delta_zeta.dat')==False):
        print('Missing spectra for: '+par+'-*step\n')
        nono+=1
    if(path.exists(filedir+suptype+'_'+par+'_plus1_'+NO+'00_cl_Delta_zeta.dat')==False):
        print('Missing spectra for: '+par+'+*step\n')
        nono+=1
    if(path.exists(filedir+suptype+'_'+par+'_minus2_'+NO+'00_cl_Delta_zeta.dat')==False):
        print('Missing spectra for: '+par+'-2*step\n')
        nono+=1
    if(path.exists(filedir+suptype+'_'+par+'_plus2_'+NO+'00_cl_Delta_zeta.dat')==False):
        print('Missing spectra for: '+par+'+2*step\n')
        nono+=1
    if nono>=1:
        print('\n')
        nono=0


# Selecting relevant data columns 
  Col = [0]
  '''
#TRACER AUTO Correlation only
for i in range(4,Nsum+3+1):
          for j in range(1-(i-3),Nsum-(i-3)+1):
              Col.append((i-1)*(Nsum+3)+i+j)
  '''
#TRACER AUTO and cross Correlation
  for i in range(1,int((Nsum)**2+1)):
    Col.append(i)

  if Nsum>0:

    Tr=loadtxt(filedir+suptype+'_'+'_fiducial'+'_'+NO+'00_cl_Delta_zeta.dat',usecols = Col)
    Tr_m = []
    Tr_p = []

    if stencil==5:
      Tr_m2 = []
      Tr_p2 = []


    for par in params:#[Keys[0:nParsLess]]: 

      Tr_m.append(loadtxt(filedir+suptype+'_'+par+'_minus1_'+NO+'00_cl_Delta_zeta.dat',usecols = Col))
      Tr_p.append(loadtxt(filedir+suptype+'_'+par+'_plus1_'+NO+'00_cl_Delta_zeta.dat',usecols = Col))

      if stencil==5:
        Tr_m2.append(loadtxt(filedir+suptype+'_'+par+'_minus2_'+NO+'00_cl_Delta_zeta.dat',usecols = Col))
        Tr_p2.append(loadtxt(filedir+suptype+'_'+par+'_plus2_'+NO+'00_cl_Delta_zeta.dat',usecols = Col))

###############################################
#Background cosmological data for lmax cut-off:
###############################################
# z-dependent lmax: In its current state, the code will use a fixed maximum ell equal to the value of 'lmax', this is because an lmax based on kmax*Chi(z) might become fairly large for higher redshifts and is not conducive for testing with the minimal example. Once you are convinced the code is correctly set-up, you can turn the z-dependent lmax cut-off on, by commenting out the line 'Lmax_Tr = [lmax for z in zmeans]' below and uncommenting the line 'Lmax_Tr = [int(kmax*Chif(z)) for z in zmeans]'. 

  Z,Ch,Hub,D,f  = np.genfromtxt(SpectraPath+suptype+'__'+'fiducial_'+NO+'00_background.dat',unpack=True,usecols = [0,4,3,17,18])

  Dofz = interpolate.interp1d(Z, D, kind='cubic',axis = 0)    #Mpc
  fofz = interpolate.interp1d(Z, f, kind='cubic',axis = 0)    #Mpc
  Chif = interpolate.interp1d(Z, Ch, kind='cubic',axis = 0)   #Mpc
  Hubf = interpolate.interp1d(Z, Hub, kind='cubic',axis = 0)  #/Mpc

  if WHICHTRACERS[0]==0 or WHICHTRACERS[0]==1:     # Selecting which surveys get lmax cut-off #(0: EUC_SPEC) #(0: SKA_SPEC)
     #Lmax_Tr = [lmax for z in zmeans]              # fixed maximum ell.
     Lmax_Tr = [int(kmax*Chif(z)) for z in zmeans]# z-dependent maximum ell.
  Lmax_Tr+=Lmax_Tr                                 # These are the lmax values for the \Delta bins then the \zeta bins (The same for each bin of corresponding redshift).
  print(Lmax_Tr)


##################
#SET SKY-FRACTION:
##################

  SKYFRAC = fskyTr#[WHICHTRACERS[0]] since without multi-tracer [TBC] there is only one experiment (and thus one sky fraction) per forecast


  NshotD=zeros([int(Nsum/2),len(N)]) 
  NshotZ=zeros([int(Nsum/2),len(N)])  
  NshotZD=zeros([int(Nsum/2),len(N)]) 

  def fi(z,zsig,zm):
    return np.exp(-1/2*(z-zm)**2/zsig**2)/np.sqrt(2*np.pi)/zsig*Chif(z)**2/Hubf(z)
  def Wi(z,zsig,zm):
    return np.exp(-1/2*(z-zm)**2/zsig**2)/np.sqrt(2*np.pi)/zsig


  def DeltaIntgd(z,zsig,zm):
     return Wi(z,zsig,zm)**2          #already normalized
# For the ranges in the normalizing integrals of the following integrand function, see the 'integral note'. 
  def ZetaIntgd(z,zsig,zm,zbar):
   if (subctype[Tracer]=='Euc_spec'):
     return (z-zbar)**2*(fi(z,zsig,zm)/quad(fi, max([0,zm-20*zsig]), zm+20*zsig,args = (zsig,zbar))[0])**2
   elif (subctype[Tracer]=='SKA_2_spec'):
     return (z-zbar)**2*(fi(z,zsig,zm)/quad(fi, max([0,zm-5*zsig]), zm+5*zsig,args = (zsig,zbar))[0])**2

  def ZetaDeltaIntgd(z,zsig,zm,zbar):
   if (subctype[Tracer]=='Euc_spec'):
     return (z-zbar)*(fi(z,zsig,zm)/quad(fi, max([0,zm-20*zsig]), zm+20*zsig,args = (zsig,zbar))[0])*Wi(z,zsig,zm)
   elif (subctype[Tracer]=='SKA_2_spec'):
     return (z-zbar)*(fi(z,zsig,zm)/quad(fi, max([0,zm-5*zsig]), zm+5*zsig,args = (zsig,zbar))[0])*Wi(z,zsig,zm)

  def ZetaNoise(zsig,zm,zbar):
   if (subctype[Tracer]=='Euc_spec'):
     return quad(ZetaIntgd, max([0,zm-20*zsig]), zm+20*zsig, args = (zsig,zm,zbar))[0]/quad(DeltaIntgd, max([0,zm-20*zsig]), zm+20*zsig,args = (zsig,zm))[0]
   elif (subctype[Tracer]=='SKA_2_spec'):
     return quad(ZetaIntgd, max([0,zm-5*zsig]), zm+5*zsig, args = (zsig,zm,zbar))[0]/quad(DeltaIntgd, max([0,zm-5*zsig]), zm+5*zsig,args = (zsig,zm))[0]


  def ZetaDeltaNoise(zsig,zm,zbar):
   if (subctype[Tracer]=='Euc_spec'):
     return quad(ZetaDeltaIntgd, max([0,zm-20*zsig]), zm+20*zsig, args = (zsig,zm,zbar))[0]/quad(DeltaIntgd, max([0,zm-20*zsig]), zm+20*zsig,args = (zsig,zm))[0]
   elif (subctype[Tracer]=='SKA_2_spec'):
     return quad(ZetaDeltaIntgd, max([0,zm-5*zsig]), zm+5*zsig, args = (zsig,zm,zbar))[0]/quad(DeltaIntgd, max([0,zm-5*zsig]), zm+5*zsig,args = (zsig,zm))[0]


  def dNdzW(z,Tracer,zsig,zm):
      return dNdz(z,Tracer)#*Wi(z,zsig,zm)

  z15 = [0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95]  # Redshifts values from 1509.07562, Table 3.
  nz15 = [6.2e-2,3.63e-2,2.16e-2,1.31e-2,8.07e-3,5.11e-3,3.27e-3,2.11e-3,1.36e-3,8.07e-4,5.56e-4,3.53e-4,2.22e-4,1.39e-4,8.55e-5,5.2e-5,3.12e-5,1.83e-5,1.05e-5] # Corresponding number density values [Mpc^-3] from 1509.07562, Table 3.
  nzf = interpolate.interp1d(z15, nz15, kind='cubic',axis = 0) #Mpc^-3

  for Tracer in WHICHTRACERS:

###############################
#NOISE AND COVARIANCE
###############################

      for i in range(len(bins)):
          zi=bins[i][0]
          zf=bins[i][1]
          zsig = zsigmas[i]
          zm = mean([zi,zf])
          if(ctype[Tracer]=='gal_spec'):
              if (subctype[Tracer]=='Euc_spec'):
                 Ngal=quad(dNdzW,zi,zf,args=(Tracer,zsig,zm))[0]
              elif (subctype[Tracer]=='SKA_2_spec'):
                 Ngal = nzf(zm)*(4*np.pi/3*(Chif(zm+zsig/2)**3-Chif(zm-zsig/2)**3))/(4*np.pi)*(np.pi/180)**2
              NshotD[i,Tracer]=1/(Ngal*(180/pi)**2)                                            # Ngal ##in /deg^2
                          
              NshotZ[i,Tracer]=1/(Ngal*(180/pi)**2)*ZetaNoise(zsigmas[i],zs[i],zbars[i])       # Ngal ##in /deg^2

              NshotZD[i,Tracer]=1/(Ngal*(180/pi)**2)*ZetaDeltaNoise(zsigmas[i],zs[i],zbars[i]) # Ngal ##in /deg^2
             


  def covDeZe(l):
# This function generates the covariance matrix for the full combination of \Delta \& \zeta at a particular ell.

      Used = []  #Since certain bins are no longer used after their cut-off lmax values, the code keeps track and only uses the necessary bins (bin indices stored in 'Used'). 
      for i in range(Nsum):
           if l<=Lmax_Tr[i] and l>=Lmin_Tr[i]:
                  Used.append(i)
      NsumTemp = len(Used)

      C=zeros([NsumTemp,NsumTemp])
    
      #Noise added to Covariance matrix (If Noise==1)
      for Tracer in WHICHTRACERS:
          NprevC = int(sum(N[0:Tracer]))
          iProxy = 0
          for i in Used:
             if i<int(sum(N[0:Tracer+1])) and i>=int(sum(N[0:Tracer])): #I.e. ONLY THE USED BINS OF THAT TRACER (bins might not be used if the lmax for that redshift has been reached.)
                if(Noise==1):
                     if iProxy<int(NsumTemp/2):
                        C[iProxy,iProxy]+=NshotD[i-NprevC,Tracer]  
                        C[iProxy,iProxy+int(NsumTemp/2)]+=NshotZD[i-NprevC,Tracer]  
                     elif iProxy>=int(NsumTemp/2):
                        C[iProxy,iProxy]+=NshotZ[int(i-Nsum/2-NprevC),Tracer]  
                        C[iProxy,iProxy-int(NsumTemp/2)]+=NshotZD[int(i-Nsum/2-NprevC),Tracer]  
                iProxy +=1 

      iProxy = 0

      for i in Used:
        jProxy = 0
        for j in Used:
             C[iProxy,jProxy]+=Tr[l-2,Nsum*(i)+(j)+1]
             jProxy +=1 
        iProxy +=1 

      return C

##########
  def covDexZe(l):
# This function generates the covariance matrix for the combination of \Delta \& \zeta at a particular ell, WITHOUT the \Delta\zeta cross correlations (dizj).
      Used = []
      for i in range(Nsum):
           if l<=Lmax_Tr[i] and l>=Lmin_Tr[i]:
                  Used.append(i)
      NsumTemp = len(Used)

      C=zeros([NsumTemp,NsumTemp])
      for Tracer in WHICHTRACERS:
          NprevC = int(sum(N[0:Tracer]))
          iProxy = 0
          for i in Used:
             if i<int(sum(N[0:Tracer+1])) and i>=int(sum(N[0:Tracer])): #I.e. ONLY THE USED BINS OF THAT TRACER (bins might not be used if the lmax for that redshift has been reached.)
                if(Noise==1):
                     if iProxy<int(NsumTemp/2):
                        C[iProxy,iProxy]+=NshotD[i-NprevC,Tracer]  
                        C[iProxy,iProxy+int(NsumTemp/2)]+=NshotZD[i-NprevC,Tracer] 
                     elif iProxy>=int(NsumTemp/2):
                        C[iProxy,iProxy]+=NshotZ[int(i-Nsum/2-NprevC),Tracer]  
                        C[iProxy,iProxy-int(NsumTemp/2)]+=NshotZD[int(i-Nsum/2-NprevC),Tracer]  
                iProxy +=1 
      iProxy = 0
      for i in Used:
        jProxy = 0
        for j in Used:
             #OFF DIAGONAL BLOCK
             if(((iProxy<NsumTemp and iProxy>=NsumTemp/2)and(jProxy<NsumTemp/2 and jProxy>=0))or(((jProxy<NsumTemp and jProxy>=NsumTemp/2)and(iProxy<NsumTemp/2 and iProxy>=0)))): 
                  C[iProxy,jProxy] = 0                                     #NO CROSS DZ, ZD
             else:
                  C[iProxy,jProxy]+=Tr[l-2,Nsum*(i)+(j)+1]#*2*pi/(l*(l+1)) #ONLY AUTO ZZ, DD
             jProxy +=1 
        iProxy +=1 
      return C

###################
#DERIVATIVE MATRIX
###################
  def deriv_covDeZe(l,p):
# This function generates the derivatives of the elements of the covariance matrix, with respect to parameter p, for the full combination of \Delta \& \zeta at a particular ell.
      
      NotUsed = []
      Used = []
      for i in range(Nsum):
           if l<=Lmax_Tr[i] and l>=Lmin_Tr[i]:
                  Used.append(i)
           else:
                  NotUsed.append(i)
      NsumTemp = len(Used)
      D=zeros([NsumTemp,NsumTemp])
      iProxy = 0
      for i in Used:
        jProxy = 0
        for j in Used:
             D[iProxy,jProxy]=derivTr(Nsum*(i)+(j)+1,p,l) 
             jProxy +=1 
        iProxy +=1 
      return D

##########
  def deriv_covDexZe(l,p):
# This function generates the derivatives of the elements of the covariance matrix, with respect to parameter p, for combination of \Delta \& \zeta at a particular ell, WITHOUT the \Delta\zeta cross correlations (dizj).
      
      NotUsed = []
      Used = []
      for i in range(Nsum):
           if l<=Lmax_Tr[i] and l>=Lmin_Tr[i]:
                  Used.append(i)
           else:
                  NotUsed.append(i)
      NsumTemp = len(Used)
      D=zeros([NsumTemp,NsumTemp])
      iProxy = 0
      for i in Used:
        jProxy = 0
        for j in Used:
             if(((iProxy<NsumTemp and iProxy>=NsumTemp/2)and(jProxy<NsumTemp/2 and jProxy>=0))or(((jProxy<NsumTemp and jProxy>=NsumTemp/2)and(iProxy<NsumTemp/2 and iProxy>=0)))): 
                  D[iProxy,jProxy] = 0						 #NO CROSS DZ, ZD 
             else:
                  D[iProxy,jProxy] = derivTr(Nsum*(i)+(j)+1,p,l)#*2*pi/(l*(l+1)) #ONLY AUTO ZZ,DD
             jProxy +=1 
        iProxy +=1 
      return D

  '''
##########################
#DERIVATIVE CHECK
##########################
# This function can be used and adapted to check the smoothness of the generated numerical derivatives, for various ells, correlations and parameters.
for p in range(len(params)):
  L = []
  Deriv = np.zeros([Nsum,Nsum,Lmax_Tr[0]-Lmin_Tr[0]])
  for ell in arange(min(Lmin_Tr),max(Lmax_Tr)):
      L.append(ell)
      NotUsed = []
      Used = []
      for i in range(Nsum):
           if ell<=Lmax_Tr[i] and ell>=Lmin_Tr[i]:
                  Used.append(i)
           else:
                  NotUsed.append(i)
      iProxy=0
      for i in Used:
        jProxy = 0
        for j in Used:
                  Deriv[iProxy,jProxy,int(ell-Lmin_Tr[0])] = derivTr(Nsum*(i)+(j)+1,p,ell)
                  jProxy +=1 
        iProxy +=1
  plt.plot(L,Deriv[0,0,:],label=Keys[p])
plt.legend()
plt.show()
  '''


#####################
#FISHER CALCULATION
#####################

  def Fkernel_DeZe(i,j,l):
    Used = []
    for k in range(Nsum):
           if l<=Lmax_Tr[k] and l>=Lmin_Tr[k]:
                  Used.append(k)
    NsumTemp = len(Used)
    ndiff = int((Nsum-NsumTemp)/2)

    covDZl = covDeZe(l)
    drvDZli = deriv_covDeZe(l,i)
    drvDZlj = deriv_covDeZe(l,j)
    FACTOR = SKYFRAC*((2*l+1)/2)
    inversa=np.linalg.inv(covDZl)
    inversaD=np.linalg.inv(covDZl[0:int(NsumTemp/2),0:int(NsumTemp/2)]) #Delta alone
    inversaZ=np.linalg.inv(covDZl[int(NsumTemp/2):int(NsumTemp),int(NsumTemp/2):int(NsumTemp)]) #Zeta alone
    drvDZliD = drvDZli[0:int(NsumTemp/2),0:int(NsumTemp/2)]
    drvDZljD = drvDZlj[0:int(NsumTemp/2),0:int(NsumTemp/2)]
    drvDZliZ = drvDZli[int(NsumTemp/2):int(NsumTemp),int(NsumTemp/2):int(NsumTemp)]
    drvDZljZ = drvDZlj[int(NsumTemp/2):int(NsumTemp),int(NsumTemp/2):int(NsumTemp)]

    inversaZD=np.linalg.inv(covDexZe(l)) #Delta and Zeta, no DeltaxZeta
    

    return FACTOR*sum(diag(mat(inversa)*mat(drvDZli)*mat(inversa)*mat(drvDZlj))), FACTOR*sum(diag(mat(inversaD)*mat(drvDZliD)*mat(inversaD)*mat(drvDZljD))), FACTOR*sum(diag(mat(inversaZ)*mat(drvDZliZ)*mat(inversaZ)*mat(drvDZljZ))), FACTOR*sum(diag(mat(inversaZD)*mat(deriv_covDexZe(l,i))*mat(inversaZD)*mat(deriv_covDexZe(l,j))))

  def F_DeZe(i,j):
    return sum([Fkernel_DeZe(i,j,l)[0] for l in arange(min(Lmin_Tr),max(Lmax_Tr)+1)]), sum([Fkernel_DeZe(i,j,l)[1] for l in arange(min(Lmin_Tr),max(Lmax_Tr)+1)]), sum([Fkernel_DeZe(i,j,l)[2] for l in arange(min(Lmin_Tr),max(Lmax_Tr)+1)]), sum([Fkernel_DeZe(i,j,l)[3] for l in arange(min(Lmin_Tr),max(Lmax_Tr)+1)])

  FisherDeZe=zeros([len(params),len(params)])
  FisherDeZeD=zeros([len(params),len(params)])
  FisherDeZeZ=zeros([len(params),len(params)])
  FisherDeZeDZ=zeros([len(params),len(params)])

  for i in range(len(params)):
    for j in range(len(params)):
        if j>=i:
            FisherDeZe[i,j]=FisherDeZe[j,i]=F_DeZe(i,j)[0]
            FisherDeZeD[i,j]=FisherDeZeD[j,i]=F_DeZe(i,j)[1]
            FisherDeZeZ[i,j]=FisherDeZeZ[j,i]=F_DeZe(i,j)[2]
            FisherDeZeDZ[i,j]=FisherDeZeDZ[j,i]=F_DeZe(i,j)[3]

  errDeZe=sqrt(diag(np.linalg.inv(FisherDeZe)))
  errDeZeD=sqrt(diag(np.linalg.inv(FisherDeZeD)))
  errDeZeZ=sqrt(diag(np.linalg.inv(FisherDeZeZ)))
  errDeZeDZ=sqrt(diag(np.linalg.inv(FisherDeZeDZ)))

  errDeZe = errDeZe[0:len(errDeZe)-Marg]
  errDeZeD = errDeZeD[0:len(errDeZeD)-Marg]
  errDeZeZ = errDeZeZ[0:len(errDeZeZ)-Marg]
  errDeZeDZ = errDeZeDZ[0:len(errDeZeDZ)-Marg]

  conderrDeZe = zeros(len(errDeZe))
  conderrDeZeD = zeros(len(errDeZeD))
  conderrDeZeZ = zeros(len(errDeZeZ))
  conderrDeZeDZ = zeros(len(errDeZeDZ))

  if Cond==1:
      for i in range(len(FisherDeZe)-Marg):
        conderrDeZe[i] = 1/sqrt(FisherDeZe[i,i])
      for i in range(len(FisherDeZeD)-Marg):
        conderrDeZeD[i] = 1/sqrt(FisherDeZeD[i,i])
      for i in range(len(FisherDeZeZ)-Marg):
        conderrDeZeZ[i] = 1/sqrt(FisherDeZeZ[i,i])
      for i in range(len(FisherDeZeDZ)-Marg):
        conderrDeZeDZ[i] = 1/sqrt(FisherDeZeDZ[i,i])

#Saving of results
  addi = ''
  if (Marg!=0):
      addi += '_marg'
  if (Fix!=0):
      addi += '_fixd'
    
  BINNES = ''
  for i in WHICHTRACERS:
      BINNES+=str(int(N[i]/2))+'bin'+subctype[i]

  RESULTS = [FisherDeZe,FisherDeZeD,FisherDeZeZ,FisherDeZeDZ]
  RESULTS2 = [errDeZe,errDeZeD,errDeZeZ,errDeZeDZ]
  RESULTS3 = [conderrDeZe,conderrDeZeD,conderrDeZeZ,conderrDeZeDZ]
  LAB2 = [suptype+addi,'Delta'+addi,'Zeta'+addi,'Delta+zeta'+addi]
  j=-1
  for LAB in ['cov','covDelta','covZeta','covDelta+zeta']:
    j+=1
    if (Noise==0):
      f=open(SpectraPath+LAB+suptype+addi+'_'+BINNES+'NN.dat','wb')
    elif (Noise==1):
      f=open(SpectraPath+LAB+suptype+addi+'_'+BINNES+'.dat','wb')
    savetxt(f,np.linalg.inv(RESULTS[j])[0:len(params)-Marg,0:len(params)-Marg])
    f.close()

    if (Noise==0):
      f=open(SpectraPath+LAB2[j]+'_fisher_'+BINNES+'NN.dat','w')
    elif (Noise==1):
      f=open(SpectraPath+LAB2[j]+'_fisher_'+BINNES+'.dat','w')
    headerline = '#Fiducial model ('
    for key in Keys[0:nParsLess]: 
        headerline+=key+','

    
    #Marginalisation component [TBC].
    if (Fix!=0):
        Marg = Fix
    if (Marg!=1):
      for i in range(2):
        if (i==2-Marg):
          headerline = headerline +','+ 'b'+str(i+1) 
        else:
          headerline = headerline+',' + 'b'+str(i+1) 
    

    f.write(headerline+'): \n')
    for i in range(len(params)-Marg):
      f.write(str(round(params[Keys[i]],10))+' ')
    
    f.write('\n \n')    
    f.write('#Marginalised errors: \n')
    for i in range(len(params)-Marg):
      f.write(str(round(RESULTS2[j][i],10))+' ')  

    f.write('\n \n')    
    f.write('#Conditional errors: \n')
    for i in range(len(params)-Marg):
      f.write(str(round(RESULTS3[j][i],10))+' ') 

    f.close()

  print('Elapsed time: ', time.time()-start, 's')
