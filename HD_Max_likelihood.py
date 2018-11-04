import struct as st
import numpy as np
import random as rd
import os as os
import scipy.io as sio
import scipy as sc


rd.seed()
#need one file per cell in it's own directory for this code to run

def doheader(fh,nlines):
	#this function will scan through the header information.
	#number of lines in header denoted by nlines
    for i in range(nlines):
        line = fh.readline()
        #getting rid of the part of the header that has the 4 byte dir listed


        if '%SAMPLING_INTERVAL' in line:
            name, value  = line.split(' ')
            samp = value[:-1]
        if 'FIRST_TIMESTAMP' in line:
            name, value  = line.split(' ')
            timestep= value[:-1]
    return float(samp), float(timestep)

    
    
#TSDpath = '/Users/stephenkeeley/Google Drive/Research/HeadDirection_2Frame/Weird/MLE_TSDs_Mar16/0930_SA1_TSD'
#DIRpath = '/Users/stephenkeeley/Google Drive/Research/HeadDirection_2Frame/Weird/MLE_TSDs_Mar16/0930_SA1_TSD'
#OUTpath = '/Users/stephenkeeley/Google Drive/Research/HeadDirection_2Frame/'

TSDpath = '/Users/Documents//TSD/SA2' #this is the folder for my testing data
DIRpath = '/Users/Documents//DIR/SA1' #this is the folder for my training data
OUTpath = '/Users/Documents//OUTPUT/' #this is my output directory

printpath = OUTpath+'fileName' 


#pull out total number of cells by looking at number of files in the directory.
#each cell must have it's own TSD file
num_cells = len(os.listdir(TSDpath))


num_samps = 20 #this is variable. is the resolution, in bins, of the estimator.

angs = np.linspace(0,350,36)
max_spks = 5*num_samps


    
#some constants  
timestamps_p_sec = 10000
estimation_window = 200 #ms


#assemble tuning probabilities
#this will generate the probability distributions for the head direction tuning
tuning_probs = list()
total_bins = np.zeros([num_cells,36])
tot_timespent = np.zeros([36])
spikes = [[] for x in xrange(num_cells)]
dirs = list()





for i in range(num_cells):
    a = np.zeros([36,max_spks+1])
    tuning_probs.append(a) 
     
#run through all cells to assemble tuning prob curves,     
for i in range(num_cells):

    listTD = os.listdir(TSDpath) 
    listDIR = os.listdir(DIRpath)
    # do your stuff
    if len(listTD) != len(listDIR):
        #print len(listTD), len(listDIR)
        raise Exception('make sure equal number of files in each directory')
    

    fileName = TSDpath +'/' + listTD[i]
    Dir_tune = DIRpath +'/' + listDIR[i]

    with open(Dir_tune, mode ='rb') as f:
        nlinbin = f.read(2)
        nlin = int(nlinbin)
        samp, timestep = doheader(f,nlin)
      
        
        timestamp_units_per_samp = timestamps_p_sec / samp #samp in samples per second

        #now lets read the data

        sim_times =[]

        binxys = f.read(7)
        while binxys != "":
            s_tot = 0
            summed_dir=0
            for k in range(num_samps):
                
                if len(binxys) != 7:
                    break
                x,y,d,s = st.unpack('=BBIB', binxys)
            
                #examine direction

                if d >360:
                    d = np.nan

                if d<=360:
                    summed_dir = summed_dir +np.exp(1j*d*(2*np.pi/360)) #generate the circular sum of all the head directions
                
                
                #generate total number of spikes in the estimator window from cell i
                s_tot= s_tot + s 

                bintimes = f.read(s*4) #read the remaining info
                binxys = f.read(7)    
                  
               
            spikes[i].append(s_tot)  #append total number of spikes for the duration of the sample time (numsamps)
            ang = np.angle(summed_dir) #find the angle for this sampe time
            if ang<0:
                ang = ang+2*np.pi

            bin_ang_ind = np.floor(36*ang/(2*np.pi)) #convert angle to index


            tuning_probs[i][bin_ang_ind][s_tot] = tuning_probs[i][bin_ang_ind][s_tot] + 1 #assemble the spikes in the correct angle indices. Do this for each numsamp chunk of time
            
            total_bins[i][bin_ang_ind] = total_bins[i][bin_ang_ind]+s_tot*1 #find the total amount of time spent at each angle, too


        f.close

    







#threshold for not having enough data to make an estimate
#tuning_probs[tuning_probs<10] = NaN

#normalize the curves

summed_tunes = np.asarray(np.sum(tuning_probs ,2)) 

#take only the first cell cause this array is redundant. (all cells from same recording have the same physical trajectory)
summed_tunes = summed_tunes[0] 

#find the weighted poisson sum(0 spikes is 0, 1 spike is 1, 2 spikes is 2, etc. number of each of these
#added up for each angle) we will divide this by the total amount of time spent at each angle
lamdas =  [i/summed_tunes for i in total_bins]


#here we regularize by finding the mean parameter for a Poisson distribution at each angle, and then filling in the probability for each number of spikes
# at each angle. That is, we calculate the the Poisson sufficient statistic (the mean) at angle, say,  30 degrees
#we used this mean parameter to determine the probability for each spike count (0,1,2,3..) at 30 degrees.  We do this for all angles
#This smoothes the results and allows us to fill in an expected number of spikes for numbers of spikes for which we had no data 
#(so we could assess this probability during testing, too).  We also found that this ML procedure acheived better estimation performance than
#simply throwing out spike count probabilities at counts and angles for which we had no data

k = range(max_spks+1)
tuning_probs = [[[np.exp(-lamdas[i][j]) * lamdas[i][j]**x / sc.math.factorial(x) for x in k]   for j in range(36)] for i in range(num_cells)]
#print tuning_probs
#tuning_probs = tuning_probs/(summed_tunes[:,None]) 

#print np.asarray(np.sum(tuning_probs,0))


#now lets estimate the position
vals2max = np.zeros(36)
dir_est = []
dir_vect = []
spk_vect= np.zeros(num_cells) 
fs = list()

#go through all the cells during test performance, read their headers
for i in range(num_cells):

    listTD = os.listdir(TSDpath) 
        
    fileName = TSDpath +'/' + listTD[i]    
    
    fs.append(open(fileName, mode = 'rb'))      

    nlinbin = fs[i].read(2)
    
    nlin = int(nlinbin)
    samp, timestep = doheader(fs[i],nlin)
    




while True:
    for j in range(num_cells):
        s_tot = 0
        for k in range(num_samps):
            binxys = fs[j].read(7)
            if len(binxys) != 7:
                break
            x,y,d,s = st.unpack('=BBIB', binxys)
            
            if j == 0:
                if d >360:
                    d = np.nan
                dir_vect.append(d)

            #if d <= 360:
            s_tot= s_tot + s
            bintimes = fs[j].read(s*4)    
            
            


        spk_vect[j] = s_tot
    #print spk_vect
    if len(binxys) != 7:
        break
    #at this range of time (numsamps) did any cells have spikes, if they did, make an estimate of the direction.
    if any(spk_vect):
        bad_est = 0
        for n in range(36):

            summ = 0
            
            for m in range(num_cells):
                
                loglik = np.log(tuning_probs[int(m)][int(n)][int(spk_vect[m])]) #find log likelihood for each position conditioned on number of spikes from a given cell

                if tuning_probs[int(m)][int(n)][int(spk_vect[m])] == 0: # should never happen now because of Poisson estimation procedure
                    loglik = -10

                summ = summ+ loglik #sum up across all cells
            vals2max[n] = summ
        
        max_ind = np.argmax(vals2max) #find index of most likely direction!
     
        ang = angs[max_ind]*(2*np.pi/360) #convert to degrees
        if ang>np.pi:
            ang = ang-2*np.pi

        dir_est.append(ang) #append the estimated angle

    else:
       dir_est.append(-10000) #if i had no spikes and couldn't make an estimate, append this number. You can throw it out in post-processing


dir_at_res = [] #we found our estimation at a fine resolution, let's now put it back to the numsamp resolution

for j in range(int(np.floor(len(dir_vect)/num_samps))):
    summed_dir=0
    for i in range(num_samps):
        if np.isnan(dir_vect[j*num_samps+i]):
                  continue
        summed_dir = summed_dir +np.exp(1j*dir_vect[j*num_samps+i]*(2*np.pi/360)) #sum it all up, take the circular mean for total number of numsamps
    if summed_dir ==0:
        dir_at_res.append(-10000)

  
    else:
        dir_at_res.append(np.angle(summed_dir))




sio.savemat(printpath, {'direction':map(float, dir_at_res),'estimate':map(float,dir_est)}) #save data



for i in range(num_cells):
    fs[i].close #close files



        
