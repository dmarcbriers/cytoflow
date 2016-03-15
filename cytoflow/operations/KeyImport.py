"""
Author: Demarcus Briers

This class allows cytoflow to automatically load the flow cytometry data, 
control experiments, and metadata of multiple flow cytometry experiments.
Note the data and metadata must be specified in a CSV file. 
An example of the CSV file,we call a KEY FILE, can be found in this direcoty.

The KEY FILE format was origionally developed by Evan Appleton.
"""

from __future__ import division
import sys
import os
import re

from collections import defaultdict
import numpy as np
import pandas as pd
import cytoflow as flow


class ImportKeyOp:
    """
    Using a STANDARD CSV file format you can easily import flow cytometry
    experiments that have different experimental conditions such as different
    media/induction concentration or different collection times. 
    
    TODO: Add support for technical replicates by looking at media/time.
  

    Attributes
    ----------
    
    keyfile : str
        path to the keyfile for the flow cytometry experiment(s)
    
    datadir : str
        path to directory/folder where ALL experimental data is located.
    
    controldir : str
        path to the directory where different control epxrimnets are.
    
    ver : str[Weiss or Cidar] DEPRECATED
        
    
    TODO
    ---------------
    Convenience Functions to implement:
    assign control dir to datadir if not set
    allow user to chose inducer
    
    """
    
    def __init__(self,keyfile,datadir="",controldir="",inducer="Dox"):
        """ Save the variables as class attributes for later."""
        #Directory information
        self.keyfile = keyfile
        self.datadir = datadir
        self.controldir = controldir
        
        #File variables
        self.experimental_data = None
        self.bleedover_controlfile = None
        self.bleedover_controlchannels = None
        
        #Save results as list then assign variables later. TODO: Make this cleaner
        self.allfiles = []
        self.allfiles = self.ImportKeyFile(keyfile,datadir,controldir,inducer)
            
        print("%i variables were saved in this object as allfiles." % len(self.allfiles))
            
    """
    def is_autofluor_control(self,row):
       
    
    def is_bead_control(self,row):
        
    
    def is_bleedover_control(self,row):
    
    
    def is_color_conversion_control(self,row):
    
    
    def is_device(self,row):

    """
    def load_negative_data(self,keyfile,datadir=""):
        """ 
        Load color control experiments used to compensate for spectral overlap
        Only returns the filename of the FIRST negative/autofluorescence control file.
        """
        #Set indices of keyfile columns
        file_col = 0
        part_col = 1
        control_col = 2
        condition_col = 3
        time_col = 4
        
        #Save data plus the driven channel for future Spectral Compensation.
        negative_experiments = []  #experiments will be added to this list
        num_rows = keyfile.shape[0]
        
        for row in range(num_rows):
            
            if  keyfile[row][control_col] == "negative":
                filename = keyfile[row][file_col]
                
                #load data
                #negative_data = flow.Tube(file = os.path.join(datadir,filename))
                negative_experiments.append(os.path.join(datadir,filename))
        
        print("Located %i autofuorescence files." % len(negative_experiments))
        return negative_experiments[0] 

    def load_color_conversion_data(self,keyfile,datadir=""):
        """ 
        Load color control experiments used to compensate for spectral overlap
        """
        #Set indices of keyfile columns
        file_col = 0
        part_col = 1
        control_col = 2
        condition_col = 3
        time_col = 4
        
        #Save data plus the driven channel for future Spectral Compensation.
        #all_experiments = []  #experiments will be added to this list
        color_conversion_files = []  #preserve order of control exp. and channel names
        num_rows = keyfile.shape[0]
        
        for row in range(num_rows):
            
            if  keyfile[row][control_col] == "color conversion":
                filename = keyfile[row][file_col]

                #load data
                #tube = flow.Tube(file = os.path.join(datadir,filename))
                color_conversion_files.append(file)
                
        ##Load data
        ### Combine files into single dataframe and import data
        #import_operation = flow.ImportOp(tubes = color_conversion_files)
        #color_control_experiments = import_operation.apply()
        
        print("Located %i color conversion files." % len(color_conversion_files))
        return color_conversion_files[0]
        #return color_control_experiments 



    def load_control_data(self,keyfile,datadir=""):
        """ 
        Functions:LLoad color control experiments used to compensate for spectral overlap.
        Input:
        Output: Dict(channel,filepath),(list) channel names
        """
        #Set indices of keyfile columns
        file_col = 0
        part_col = 1
        control_col = 2
        condition_col = 3
        time_col = 4
        
        #Save data plus the driven channel for future Spectral Compensation.
        color_controls = defaultdict(lambda:0)  #experiments will be added to this list
        control_channels = []  #preserve order of control exp. and channel names
       
        num_rows = keyfile.shape[0]
        for row in range(num_rows):
            if  keyfile[row][control_col] == "color control":
                #Parse Media conditions
                filename = keyfile[row][file_col]
                driven_channel = keyfile[row][part_col]

                #Save Color Control channels
                color_controls[driven_channel] = os.path.join(datadir,filename)
                control_channels.append(driven_channel)
        
        print("Located %i control files." % len(control_channels))
        return color_controls,control_channels 

    def load_bead_data(self,keyfile,datadir=""):
        """ 
        For bioCPS gene circuits experiments vary the inducer concentration.
        """
        #Set indices of keyfile columns
        file_col = 0
        part_col = 1
        control_col = 2
        condition_col = 3
        time_col = 4
        
        beads_experiments = [] #experiments will be added to this list
        num_rows = keyfile.shape[0]
        
        for row in range(num_rows):
            if  keyfile[row][part_col] == "beads" or keyfile[row][control_col] == "beads":
                filename = keyfile[row][file_col]
                #tube = flow.Tube(file = os.path.join(datadir,filename))
                beads_experiments.append(os.path.join(datadir,filename))
                
        ##Load data
        ### Combine files into single dataframe and import data
        #import_operation = flow.ImportOp(tubes = all_experiments)
        #bead_experiments = import_operation.apply()
        
        print("Located %i bead files." % len(beads_experiments))
        return beads_experiments[0]     


    def load_experimental_data(self,keyfile,datadir="",inducer="Inducer Level"):
        """ 
        We detect experimental file by the control column being empty
        """
        #Set indices of keyfile columns
        file_col = 0
        part_col = 1
        control_col = 2
        condition_col = 3
        time_col = 4
        regulation_col = 5
        
        all_experiments = [] #experiments will be added to this list
        #Keep track of replicates with 2-level dict with default count of 0.
        # Dict[media_conditions][time] = 0
        replicates_dict = defaultdict(lambda:defaultdict(lambda:0)) 
        
        num_rows = keyfile.shape[0]
        for row in range(num_rows):
            if  keyfile[row][control_col] == "":
                filename = keyfile[row][file_col]
                experiment_time = keyfile[row][time_col]

                #Parse Media conditions which are numbers. Ignore units like  nM or ul/M
                condition = keyfile[row][condition_col]
                condition = re.search(r'[\d|.]{1,}',condition).group()
                
                #Check for replicate experiments
                replicates_dict[str(condition)][str(experiment_time)] +=1
                rep = replicates_dict[str(condition)][str(experiment_time)]
                
                #Set the metadata for the experiment
                #print("Loading replicate: %i" % rep)
                tube = flow.Tube(file = os.path.join(datadir,filename),
                                 conditions = {inducer : float(condition),
                                               'Replicate': int(rep) }
                                )
                
                all_experiments.append(tube)
                
        ##Load data
        ### Combine files into single dataframe and import data
        import_op = flow.ImportOp(conditions = {'Dox' : "float",
                                                'Replicate' : "int"},
                                                 tubes = all_experiments,
                                                 coarse = False,
                                                 coarse_events = 10000
                                  )
        cytometry_experiments = import_op.apply()
        
        print("Successfully loaded %i experimental files." % len(all_experiments))
        return cytometry_experiments       

    def load_all_data(self,keyfile):
        """ 
        Only loop through file one time.
        """
        
        num_rows = keyfile.shape[0]
        
        all_experiments = []  #Load experimental data
        blank_controls = []   #hold autofluroescence control
        color_compensate_controls = []  # hold compensation/bleedover controls
        color_conversion_controls = []  #hold data for converting to single color
        beads_controls = []
        
        
    def ImportKeyFile(self,keyfile,datadir="",controldir="",inducer="Dox"):
        """
        Read keyfile from a STANDARD format. Will give link with an example.
        
        Keyfile
        -----------------
        Filename - absolute or relative name of FCS file.
        Part - device/level0/level1/lvel2. These are MoClo names.
        Control - beads|beads/negative|blank/color controls/color conversion.
        Condition - Media conditions of the experiment. Usually this is concentration of the inducer.
        Time - How many days/hours after transfection/transformation was data collected.
        ------------------
        
        This function also detects replicates by the same time media condition.
        """
        
        #1. Load keyfile as numpy array
        keyfile = np.genfromtxt(keyfile,
                                delimiter=',',
                                dtype=str,
                                skip_header=True)
        
        #Is this a valid keyfile with rewquired experiments
        #print("Checking if key file is valid")
        
        
        #Parse Keyfile
        print("Loading experimental data from %s" % datadir)
        experimental_data = self.load_experimental_data(keyfile,datadir)
        control_data,control_channels = self.load_control_data(keyfile,controldir)
        negative_data = self.load_negative_data(keyfile,controldir)
        bead_file = self.load_bead_data(keyfile,controldir)
        color_conversion = self.load_color_conversion_data(keyfile,controldir)
        
        return experimental_data,control_data,control_channels,negative_data,bead_file,color_conversion
