import numpy as np
import os
import csv
import subprocess
import time
import shutil
import psutil
from itertools import combinations

# generate single replications
def generateRep(rep_dir, anatomy, vowels, nGenerations, fitness, nHidden, 
                mutationRate, crossoverRate, activation, n_formants, taufactor, sigmaScaling, 
                plusSelection, rankingSelection, parent_selection, 
                offspring_selection, popSize, nElites, config_root, wav):
    
    if os.path.exists(rep_dir):
        shutil.rmtree(rep_dir)
    
    # java parameters
    j_params = (("problem","vtl"),("type","janssen"),("fitness",fitness),("activation",activation),
                ("parentSelection",parent_selection),("offspringSelection",offspring_selection),
                ("rankingSelection",rankingSelection),("plusSelection",plusSelection),
                ("sigmaScaling", sigmaScaling),("wav",wav),
                ("mseExponent",0.5), ("nIterations",nGenerations),("popSize",popSize),
                ("mutationRate",mutationRate),("crossoverRate",crossoverRate),("nThreads",1),
                ("iAnatomy",anatomy),("nFormants",n_formants),("nElites",nElites),
                ("tauFactor",taufactor),("lambdaFactor",1))
    
    
    print "\tgenerating " + rep_dir  
    
    # write dirs
    os.makedirs(os.path.join(rep_dir))

    # write java config file
    with open(os.path.join(rep_dir,"config.csv"), 'wb') as config_file:
        writer = csv.writer(config_file)
        
        for param in j_params[:10]:
            writer.writerow(param)
            
        writer.writerow(["nHidden"] + list(nHidden))
        writer.writerow(["targets"] + list(vowels))
            
        for param in j_params[10:]:
            writer.writerow(param)
            
    # copy binaries
    try:            
        for binary in os.listdir(config_root):
            source = os.path.join(config_root,binary)
            
            if "anatomy" in binary: 
                destination = os.path.join(rep_dir,"anatomy.csv")
                shutil.copy(source, destination)
            else:
                destination = os.path.join(rep_dir,binary)
                shutil.copy(source, destination)
    except OSError:
        pass

# scan directory for previously completed reps
def scan(super_path, max_generations, n_replications):
    (running,finished,latest_rep_gens) = ([],[],[])
    n_completed = 0
        
    for set in os.listdir(super_path):
        (set_finished,set_completed) = ([],[])
        
        if set[0] != "_":
            set_path = os.path.join(super_path,set)
            
            for rep in os.listdir(set_path):
                if rep[0] != "_":
                    rep_path = os.path.join(set_path,rep)
                    log_path = os.path.join(rep_path, "output.txt")
                    
                    # first check for finished marker
                    try:
                        with open(log_path, 'rb') as log_file:
                            for line in log_file:
                                pass

                            isFinished = True if line.startswith("Finished!") else False
                            isRunning = not isFinished
                            
                    # when output.txt does not exist
                    except IOError: 
                        (isFinished,isRunning) = (False,False)
                    # when output.txt is empty
                    except UnboundLocalError:
                        (isFinished,isRunning) = (False,True)
                        
                    if isFinished:
                        finished.append((set,"",rep))
                        set_finished.append((set,"",rep))
                    if isRunning:
                        running.append((set,"",rep))
            
            # determine latest reps                      
            try:       
                for rep in os.listdir(os.path.join(set_path,"_completed")):
                    set_completed.append((set,"_completed",rep))
                    n_completed += 1     
            except WindowsError:
                pass
    
            for i_rep in xrange(n_replications):
                reps = filter(lambda x: int(x[-1][3:x[-1].index(".")]) == i_rep, set_completed + set_finished)            
                reps = sorted(reps, key=lambda x: int(x[-1][x[-1].index(".")+1:]))
                
                try:
                    latest_rep_gens.append(reps[-1])
                except IndexError:
                    pass                
                   
    return (running,latest_rep_gens,n_completed)

# try to create new generation for running replications
def attemptNewGens(super_path, finisheds, n_chain_gen, n_formants, n_vowels): 
    new_gens = []
       
    for (set_dir,inter_dir,rep_dir) in finisheds:
            
        rep_src_path = os.path.join(super_path,set_dir,inter_dir,rep_dir)

        try:
            # first try if file lock is released
            output_name = os.path.join(rep_src_path, "output.txt")
            os.rename(output_name, output_name + "_temp")
            os.rename(output_name + "_temp", output_name)
                    
            # create next-gen dir
            (rep,chain_gen) = rep_dir.split(".")
            chain_gen = int(chain_gen)
    
            if chain_gen <= n_chain_gen:         
                if chain_gen < n_chain_gen:
                    next_gen = chain_gen + 1
                    follow_up_dst = os.path.join(super_path,set_dir,rep + "." + str(next_gen))
                    
                    # this also triggers the try statement
                    os.makedirs(follow_up_dst)
                    print " \t" + set_dir + ": " + rep_dir + "->" + rep + "." + str(next_gen)
                    
                    # read next-gen targets                
                    with open(os.path.join(rep_src_path, "logElitesPhenotypes.csv"), 'rb') as pheno_file:
                        pheno_reader = csv.reader(pheno_file)
                        
                        for (i_line,line) in enumerate(pheno_reader):
                            if i_line == 0:
                                header = line[n_formants*n_vowels+1:-4*n_vowels-7]
                                n_params = len(header) / n_vowels
                                header = [header[i*n_params:(i+1)*n_params] for i in xrange(n_vowels)]
                                param_names = [v[v.index("_")+1:] for v in header[0]]
                                vowel_names = [v[0][:v[0].index("_")] for v in header]
                                
                        params = line[n_formants*n_vowels+1:-4*n_vowels-1]
                        params = [params[i*n_params:(i+1)*n_params] for i in xrange(n_vowels)]
                    
                    # read next-gen anatomy
                    if next_gen == 1:
                        with open(os.path.join(rep_src_path, "config.csv"), 'rb') as config_file:
                            config_reader = csv.reader(config_file)
                            
                            for line in config_reader:
                                if line[0] == "iAnatomy":
                                    iAnatomy = int(line[1])
                                    break
                        
                        with open(os.path.join(rep_src_path, "anatomy.csv"), 'rb') as anatomy_file:
                            anatomy_reader = csv.reader(anatomy_file)
                            
                            for (i_line,line) in enumerate(anatomy_reader):
                                if i_line == iAnatomy + 2:
                                    anatomy = line[1:]
                    
                    # generate next-gen target file
                    targets_src_file = os.path.join(rep_src_path,"targets.csv")
                    targets_dst_csv = os.path.join(follow_up_dst,"targets.csv")
                    
                    with open(targets_src_file, 'rb') as src_file, open(targets_dst_csv, 'wb') as dst_file:
                        targets_reader = csv.reader(src_file)
                        targets_writer = csv.writer(dst_file)
                        
                        named_targets = zip(vowel_names,params)
                        
                        for (i_line,line) in enumerate(targets_reader):
                            if next_gen == 1 and i_line == 1:
                                targets_writer.writerow(anatomy)
                                             
                            elif line[0] in vowel_names:
                                try:
                                    i_target = vowel_names.index(line[0])
                                    (target_name,target_params) = named_targets[i_target]
                                    targets_writer.writerow([target_name] + target_params)
                                except ValueError:
                                    pass                        
                            else:
                                targets_writer.writerow(line)
                    
                    # and queue next gen
                    rep = [set_dir, rep[3:], next_gen]
                    new_gens.append(rep)
                
                    # copy over last config files
                    for item in ["config.csv","anatomy.csv"]:
                        shutil.copyfile(os.path.join(rep_src_path,item), os.path.join(follow_up_dst,item))
                
                # move finished rep to complete
                rep_dst_path = os.path.join(super_path,set_dir, "_completed", rep_dir)
                
                if rep_src_path != rep_dst_path: 
                    shutil.move(rep_src_path, rep_dst_path)
                    
        except WindowsError, e:
            pass
        
    return new_gens    

# read global config file
def readConfig():
    parameters = {}
    
    with open("config.csv", 'rb') as config_file:
        reader = csv.reader(config_file)
        
        for line in reader:
            key = line[0]
            string_values = line[1:]
            
            try:
                string_values = string_values[:string_values.index("")]
            except ValueError:
                pass
            
            values = []
            for value in string_values:
                try:
                    try:
                        value = int(value)
                    except ValueError: # if float
                        value = float(value)
                    value = int(value) if int(value) == float(value) else float(value)
                except ValueError: # if string
                    pass
                
                values.append(value)
                
            if not key=="nHidden" and len(values) == 1:
                values = values[0]                
        
            parameters[key] = values
            
    return parameters            

# main 
if __name__ == "__main__":
    # set priority to low
    proc = psutil.Process(os.getpid())
    proc.nice(psutil.IDLE_PRIORITY_CLASS)
    prio = "low"
    
    # read and set parameters
    parameters = readConfig()
    
    config_root = parameters["config_root"]
    nGenerations = parameters["nIterations"]
    n_formants = parameters["nFormants"]
    super_path = os.path.join(parameters["data_root"], str(parameters["expLabel"]))

    base_rep = [nGenerations, parameters["fitness"], parameters["nHidden"], 
    parameters["mutationRate"], parameters["crossoverRate"], parameters["activation"],
    n_formants, parameters["tauFactor"], parameters["sigmaScaling"], 
    parameters["plusSelection"], parameters["rankingSelection"], parameters["parentSelection"], 
    parameters["offspringSelection"], parameters["popSize"], parameters["nElites"], config_root,
    parameters["wav"]]
    
    # set vowels
    #'i','I','y','Y','e','E','oe','OE','ae','u','U','o','O','a','i-','schwa','er','r'
    n_targets = parameters["nTargets"]
    vowels = [list(v) for v in (combinations(parameters["targets"], n_targets))]
    
    # read anatomy file
    with open(os.path.join(config_root, "anatomy.csv")) as target_file:
        reader = csv.reader(target_file)
        anatomies = [line[0].replace(" ", "_") for line in reader][2:]
    
    # init
    n_replications = parameters["nReplications"]
    i_anatomies = parameters["iAnatomies"]
    n_chain_gen = parameters["nChainGen"]
    
    # if new run
    if(not os.path.exists(super_path)):
        os.makedirs(super_path)
        n_completed = 0
    # if continuing from previous run  
    else:             
        # remove unfinisheds, transfer finished earliers
        (unfinisheds,latest_rep_gens,n_completed) = scan(super_path, nGenerations, n_replications)
        
        for unfinished in unfinisheds:
            shutil.rmtree(os.path.join(super_path, *unfinished))   
                   
    # generate initials sets/reps (but exclude already finished before)
    reps = []   
    
    for i_anatomy in i_anatomies:
        for vowel_subset in vowels:
            set_dir = anatomies[i_anatomy] + "."
            
            for vowel in vowel_subset:
                set_dir += vowel + "_"
                
            set_dir = set_dir[:-1]
            
            for i_rep in xrange(n_replications):                       
                # only add if not computed before   
                already_computed = False
                
                try:
                    for item in latest_rep_gens:
                        (f_set_dir,f_i_rep) = (item[0],item[-1]) 
                        
                        if f_set_dir == set_dir and int(f_i_rep.split(".")[0][3:]) == i_rep:
                            already_computed = True
                            break
                # if first init
                except NameError:
                    pass
                
                if not already_computed:
                    rep = [set_dir, i_rep, 0, i_anatomy, vowel_subset] + base_rep
                    reps.append(rep)             
      
    # keep try running new reps
    n_reps = len(vowels) * len(i_anatomies) * n_replications
    
    while n_completed < n_reps * n_chain_gen:               
        (runnings,latest_rep_gens,n_completed) = scan(super_path, nGenerations, n_replications)
        
        new_reps = attemptNewGens(super_path, latest_rep_gens, n_chain_gen-1, n_formants, n_targets)     
        
        reps += new_reps
        reps = sorted(reps, key=lambda x: (-x[2], int(x[1])))
            
        # run reps if slots are available
        if len(reps) and len(runnings) < parameters["maxProcesses"]:            
            rep = reps[0]
            set_dir = str(rep[0])
            rep_dir = "rep" + str(rep[1]) + "." + str(rep[2])
            
            set_path = os.path.join(super_path, set_dir)
            rep_path = os.path.join(set_path, rep_dir)
            
            if len(rep) > 3:                
                generateRep(rep_path, *rep[3:])
             
            # submit command            
            command = "start /" + prio + " /b java -jar " + "Agent.jar " + rep_path + "\ >> " + os.path.join(rep_path,"output.txt")
            print "\trunning " + set_dir + ", " + rep_dir + " ..."
            subprocess.Popen(command, env={'PATH': parameters["java_path"]}, shell=True)
            
            # if resuming, remove preceding gen from queue 
            reps = reps[1:]            
            time.sleep(1)
        else:
            print "waiting..."
            time.sleep(60)