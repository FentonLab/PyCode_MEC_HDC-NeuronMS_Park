
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:27:09 2015

@author: stephenkeeley
"""

import struct as st
import numpy as np
import random as rd
import os as os
import scipy.io as sio
rd.seed()



    
    

def doheader(fh,nlines):
    for i in range(nlines):
        line = fh.readline()
        #getting rid of the part of the header that has the 4 byte dir listed
        #this will read the file up through the header
        if '%SAMPLING_INTERVAL' in line:
            name, value  = line.split(' ')
            samp = value[:-1]
        if 'FIRST_TIMESTAMP' in line:
            name, value  = line.split(' ')
            timestep= value[:-1]
    return float(samp), float(timestep)

    
#roots for the data, change these as necessary    
TSDpath_root = '/Users/HeadDirection_2Frame/spk_trig_15min_2/'
OUTpath = '/Users/HeadDirection_2Frame/spk_trig_5s_rate/'






for folder in os.listdir(TSDpath_root):
    if folder.endswith(".DS_Store"):
        continue

    titles  = []
    ends = []

    TSDpath = TSDpath_root +'/' + folder
    
    #where to print the output
    printpath = OUTpath+folder+'spike_triggered_tunes_5s_RATE_2nd15'

    
    
    num_cells = len(os.listdir(TSDpath))
    #num_samps = 30 #time after trigger to make a curve
    frs = np.zeros([num_cells,37])
    angs = np.linspace(0,350,36)

    
    
        
      
    timestamps_p_sec = 10000
    estimation_window = 200 #ms
    
    
    #assemble tuning probabilities
    tuning_spks =np.zeros([36,num_cells]) 
    tuning_times = np.zeros([36,num_cells])
    tuning_curves = list()
    names = list()
    
    spikes = [[] for x in xrange(num_cells)]
    angs = [[] for x in xrange(num_cells)]
    angs_all = [[] for x in xrange(num_cells)]
    times = [[] for x in xrange(num_cells)]
    times_all=[[] for x in xrange(num_cells)]

    
    threshold = 0
        

    #cycle through total cells    
    for i in range(num_cells):
    
        listTD = os.listdir(TSDpath) 
        fileName = TSDpath +'/' + listTD[i]
        titles.append(listTD[i])


        ends.append(fileName[fileName.find('P'):])
        
        with open(fileName, mode ='rb') as f:
            nlinbin = f.read(2)
            print fileName, nlinbin
    
            nlin = int(nlinbin)
            samp, timestep = doheader(f,nlin)
          
            
            
            timestamp_units_per_samp = timestamps_p_sec / samp #samp in samples per second
            
    
            #now lets read the data
    
            sim_times =[]
    
            binxys = f.read(7)
            timer = timestep
            while binxys != "":

                s_tot = 0
                summed_dir=0
    
                    
                if len(binxys) != 7:
                    break
                x,y,d,s = st.unpack('=BBIB', binxys)
                bintimes = f.read(s*4)

                if s>0:
                    #make string
                    stringg = ''
                    for x in xrange(s):
                        stringg = stringg+'I'
                    time = st.unpack(stringg, bintimes)

                #examine direction
    
                if d<=360 and s>0:
                    spikes[i].append(s)
                    angs[i].append(2*d*np.pi/360)
                    times[i].append(time[0])
                if d<=360:
                    angs_all[i].append(2*d*np.pi/360)
                    times_all[i].append(timer)
                    
                timer = timer+ timestamp_units_per_samp 

                
                binxys = f.read(7)        
            #normalize the curves
    
            f.close
            
            #now calc spike triggered ang diffs
            
    
    diff_tunes = np.zeros([num_cells,num_cells,37])#tuning curves for pairs of cells in spikes
    times_tunes = np.zeros([num_cells,37])#tuning curves for pairs of cells in timebin counts
    rate_tunes = np.zeros([num_cells,num_cells,37])#tuning curves for pairs of cells in rate

    cell_pair_info = np.zeros((num_cells*num_cells,), dtype=np.object)
    
    
    
    
    for i in range(num_cells):
        for j in range(i,num_cells):
            if i!=j:
                #generate titles
                cell_pair_info[i*num_cells+j] = titles[i]+'_'+ends[j]
    
    
    
    
    
    move_ts =900#number of timesteps to move ahead and behind on a given spike
    #diff_counter = -1 #initialize the pair of cells counter
    for i in range(num_cells):
        indexer =[0]*num_cells
        time_indexer =[0]*num_cells
        for j in xrange(len(spikes[i])):
            #go through cell and trigger spikes
            current_ang = angs[i][j]#ang at spike j in cell i
            start_time= times[i][j] - move_ts*timestamp_units_per_samp#time of that spike minus the window

            end_time = start_time + move_ts*timestamp_units_per_samp#time to move to
            
            
            
            #assemble time vect around each spike.
            times_cell_i = times_all[i]            
            while times_cell_i[time_indexer[i]]<start_time:
                
                time_indexer[i]=time_indexer[i]+1
                
            indexgt_trig = time_indexer[i]
                
            while times_cell_i[indexgt_trig]<end_time:
                #bin Times here!!!
            
                ang_diff = current_ang- angs_all[i][time_indexer[i]]
                #print ang_diff#, i,k,indexer,timevect[indexgt_trig],end_time
                if ang_diff >=np.pi:
                    ang_diff = ang_diff-2*np.pi
                if ang_diff <-np.pi:
                    ang_diff = ang_diff+2*np.pi  
            
            
            
                bin_ang_ind = np.floor(36*(ang_diff+np.pi)/(2*np.pi))
                times_tunes[i][bin_ang_ind] = times_tunes[i][bin_ang_ind] +1
                
                
                indexgt_trig = indexgt_trig+1
                if indexgt_trig>=len(times_cell_i):
                    break
                time_indexer[i]=time_indexer[i]+1

                
                
            for k in range(num_cells):
                if k!=i:
                    #go through all other cells
                    #find times where n steps into the future these other guys spiked
                    timevect = times[k] #spiketimes for other cells spikes
                    #looking for the times that fall between start time and end time
                    
                    if indexer[k]<len(timevect):
                        while timevect[indexer[k]]<start_time:
                            #break when end of vect for max value here
                            indexer[k]=indexer[k]+1
                            if indexer[k]>=len(timevect):
                                break
                        indexgt_trig = indexer[k]

                        if indexgt_trig>=len(timevect):
                            break
            
                        while timevect[indexgt_trig]<end_time:
                            #break when end of vect for max value here
            
            
                            #record all the stuff
                            ang_diff = current_ang- angs[k][indexgt_trig]

                            #keep angle between neg pi and pi
                            if ang_diff >=np.pi:
                                ang_diff = ang_diff-2*np.pi
                            if ang_diff <-np.pi:
                                ang_diff = ang_diff+2*np.pi    
                                

                            bin_ang_ind = np.floor(36*(ang_diff+np.pi)/(2*np.pi))
                            diff_tunes[i][k][bin_ang_ind] = diff_tunes[i][k][bin_ang_ind] +1

                            indexgt_trig = indexgt_trig+1
                            if indexgt_trig>=len(timevect):
                                break
                            
                            
            rate_tunes[i] = diff_tunes[i]/times_tunes[i][None,:]
                            
                            
                            


    #save output
    sio.savemat(printpath, {'diff_tunes':diff_tunes,'rate_tunes':rate_tunes, 'times_tunes':times_tunes,'cell_info':cell_pair_info})
                                


