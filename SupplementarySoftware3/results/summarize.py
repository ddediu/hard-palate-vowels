import os
import csv
import xml.etree.ElementTree as ET
from scipy.interpolate import interp1d

def getValues(line, n_targets, n_formants, n_global_nishimuras, n_var_nishimuras):
    i_params = n_formants*n_targets+1
        
    if line[-1] == '':
        line = line[:-1]
        
    try:
        line[-1].index("_")
        n_nishimuras = 0
    except ValueError:
        n_nishimuras = n_global_nishimuras + n_var_nishimuras * n_targets
                             
    generation = line[0]
    
    formants = line[1:i_params]
    formants = [formants[i*n_formants:(i+1)*n_formants] for i in xrange(n_targets)]
    
    if n_nishimuras == 0:
        params = line[i_params:]
    else:
        params = line[i_params:-n_nishimuras]
        
    n_params = len(params) / n_targets
    params = [params[i*n_params:(i+1)*n_params] for i in xrange(n_targets)]
    
    nishimuras = line[-n_nishimuras:]
    nishimuras = [nishimuras[i*n_var_nishimuras:(i+1)*n_var_nishimuras] for i in xrange(n_targets)]
    
    if n_global_nishimuras != 0:
        nishimuras = [v + line[-n_global_nishimuras:] for v in nishimuras]
       
    return (generation,formants,params,nishimuras)   

def getFixedParam(rep_path, i_anatomy):
    abbreviations = {"Vert. hyoid pos.":"HY", 'Horz. jaw pos.':"JX", 'Velum shape':"VS", 
                     'Velic opening':"VO", 'Wall compliance':"WC", 'Tongue side elevation 1':"TS1", 
                     'Tongue side elevation 2':"TS2", 'Tongue side elevation 3':"TS3", 
                     'Tongue side elevation 4':"TS4", 'Min. area (tongue body)':"MA1", 
                     'Min. area (tongue tip)':"MA2", 'Min. area (teeth-lips)':"MA3", 
                     'Maxilla width':"HPZ", 'Maxilla length':"HPX", 'Maxilla curvature':"HPC", 
                     'Mandible width':"JAZ", 'Mandible length':"JAX", 'Mandible curvature':"JAC", 
                     'Maxilla angle':"HPA", 'Palate fronting':"PA1", 'Palatal concavity':"PA2", 
                     'Palate angle':"PA3", 'Alveo-palatal weight':"PA4", "SVTv length":"LEN"}
    
    with open(os.path.join(rep_path,"anatomy.csv")) as anatomy_file:
        reader = csv.reader(anatomy_file)
        
        reader.next()        
        header_abbrevs = [abbreviations[item] for item in reader.next()[1:]] 
        
        for (i_line,line) in enumerate(reader):
            if i_line == i_anatomy:
                break
                        
    return (header_abbrevs,line[1:])

def terminationCondition(timeseries, n_generations, ax1=None, ax2=None):
    window_size = n_generations / 5
    
    derivatives = []
    
    for i in xrange(len(timeseries) - window_size):
        dy = timeseries[i+window_size] - timeseries[i]
        derivative = dy / window_size
        derivatives.append(derivative)
    
    xs = range(window_size,len(timeseries))       
    
    threshold = 0
    belows = [i for (i,v) in enumerate(derivatives) if v >= threshold]
    
    try:
        ax2.plot(xs, derivatives, label="window: " + str(window_size), lw=1)
        ax2.axhline(y=threshold, color='k', lw=1, ls='--')
    except AttributeError:
        pass
    
    try:        
        terminator_x = belows[0]
        terminator_y = timeseries[terminator_x]
        
        try:
            ax1.axvline(x=terminator_x, color='k', lw=1, ls='--')
            ax1.axvline(x=terminator_x + window_size, color='k', lw=1, ls='--')
            ax2.axvline(x=terminator_x + window_size, color='k', lw=1, ls='--')
        except AttributeError:
            pass
            
        terminator = (terminator_x,terminator_y)  
    except IndexError:
        terminator = None
    
    return terminator

def get_rep_data(repRoot):    
    with open(os.path.join(repRoot,"logPopulation.csv"),'rb') as populationFile:
        nGenerations = len(list(populationFile)) - 1
    
    elite_data = []
    with open(os.path.join(repRoot,"logElitesGenotypes.csv"),'rb') as elitesFile:        
        eliteReader = csv.reader(elitesFile)
        eliteHeader = eliteReader.next()
        
        for (i,field) in enumerate(eliteHeader):
            if field.startswith("L0"):
                error_titles = eliteHeader[1:i]
                i_weights_start = i
                break
            
        for (i,field) in enumerate(eliteHeader):            
            if field.startswith("step"):
                i_step_start = i
                break
              
        for line in eliteReader:
            values = [float(v) for v in line[:-1]]
            
            if values[0] >= 0:
                elite_data.append(values)
                
    elite_data = zip(*elite_data)
    (generations,errors) = (elite_data[0],elite_data[1:i_weights_start])
    
    xs = [x for x in xrange(int(generations[-1])+1)]               
    
    interpolated_errors = []
    for error in errors:
        interpolator = interp1d(generations, error, kind='zero')
        interpolated_error = interpolator(xs)
        
        interpolated_error = list(interpolated_error)    
        interpolated_error += [interpolated_error[-1]] * ((nGenerations + 1) - len(interpolated_error))        
        
        interpolated_errors.append(interpolated_error)       
    
    return (interpolated_errors,nGenerations)


def getElite(set_path,rep,n_maxed_out,use_x_terminator=True):
    rep_path = os.path.join(set_path,rep)
    
    with open(os.path.join(rep_path,"config.csv"), 'rb') as config_file:
        reader = csv.reader(config_file)
        
        for line in reader:
            if line[0] == "targets":
                vowels = line[1:]
                n_targets = len(vowels)
            elif line[0] == "nFormants":
                n_formants = int(line[1])
            elif line[0] == "iAnatomy":
                i_anatomy = int(line[1])
      
    (fixed_header,fixed_param) = getFixedParam(rep_path, i_anatomy)
    (interpolated_errors,max_generation) = get_rep_data(rep_path)
    
    if use_x_terminator:
        try:
            (x_terminator,_) = terminationCondition(interpolated_errors[0], max_generation+1)
        except TypeError:
            x_terminator = max_generation + 1
            n_maxed_out += 1
    else:
        x_terminator = max_generation + 1
    
    elite_path = os.path.join(rep_path,"logElitesPhenotypes.csv")
    
    with open(elite_path, 'rb') as elite_file:
        reader = csv.reader(elite_file)
        
        headers = [reader.next() for _ in xrange(3)]
        headers = [item[:-1] if item[-1]== '' else item for item in headers]
        (header,target,alt) = headers
        
        for line in reader:
            if line[-1] == '':
                line = line[:-1]
                                      
            try:
                generation = int(float(line[0]))
            
                if generation > x_terminator:
                    break
                else:
                    n_line = line
            except ValueError:
                pass
    
    if len(header) > len(alt):
        alt = [alt[0]] + [0] * (len(header) - len(alt)) + alt[1:]
        
    return (fixed_header,fixed_param,header,target,alt,n_line,n_targets,n_formants,vowels,n_maxed_out)

def writeXML(base_tree, vowel, header_fixed, params_fixed, header_params, values_params):
    root_node = base_tree.getroot()
    node_vt = root_node.findall("vocal_tract_model")[0]
    node_shapes = node_vt.findall("shapes")[0]
    
    new_shape = ET.Element("shape")
    new_shape.set("name", vowel)
    
    header = header_fixed + header_params
    values = params_fixed + values_params
    
    for (name,param) in zip(header,values):
        new_param = ET.Element("param")
        new_param.set("value", param)
        new_param.set("name", name)
        new_param.set("domi", "100.0")
        
        new_shape.append(new_param)
    
    node_shapes.append(new_shape)

if __name__ == "__main__":
    root = r"./"
    
    #i_truncate = 50
    i_truncate = None
    use_x_terminator = False
    
    sets = os.listdir(root)
    output = []
    
    tree = ET.parse("JD2.speaker")
    n_maxed_out = [0] * len(sets)
    
    for (i_set,set) in enumerate(sets):
        if set[0] != "_":
            try:
                set_path = os.path.join(root, set, "_completed")
                reps = os.listdir(set_path)
            except WindowsError:
                reps = []
            
            reps = sorted(reps, key=lambda x: int(x[3:x.index(".")]))      
            
            try:
                reps = reps[:i_truncate*50]
            except TypeError:
                pass
            
            for (irep,rep) in enumerate(reps):
                if rep[0] != "_":  
                    # get elite data
                    (header_fixed,params_fixed,header,target,alt,elite,n_targets,n_formants,vowels,n_maxed_out[i_set]) = getElite(set_path,rep,n_maxed_out[i_set],use_x_terminator)
                    
                    n_global_nishimuras=6
                    n_var_nishimuras=4
                    
                    # get values and generate header
                    (_,header_formants,header_params,header_nishimura) = getValues(header, n_targets, n_formants, n_global_nishimuras, n_var_nishimuras)
                    
                    new_header_nishimura = []
                    for vowel in header_nishimura:                    
                        try:
                            vowel = [v[v.index("_")+1:] for v in vowel[:n_var_nishimuras]] + vowel[-n_global_nishimuras:]
                        except ValueError:
                            pass
                        
                        new_header_nishimura.append(vowel) 
                    header_nishimura = new_header_nishimura
                    
                    if len(output) == 0:
                        output_header = ["condition","vowel","replication","chain_gen"]
                        
                        for data_type in ("_target","_alt","_elite"):
                            output_header += [item[item.index("_")+1:] + data_type for item in header_formants[0]]
                            output_header += [item[item.index("_")+1:] + data_type for item in header_params[0]]
                            output_header += [item + data_type for item in header_nishimura[0]]
                        
                        output_header.append("generation")
                        
                        output.append(output_header)           
                    
                    # get targets, alts and elite                 
                    (_,targets_formants,targets_params,targets_nishimura) = getValues(target, n_targets, n_formants, n_global_nishimuras, n_var_nishimuras)                    
                    (_,alts_formants,alts_params,alts_nishimura) = getValues(alt, n_targets, n_formants,  n_global_nishimuras, n_var_nishimuras)   
                    (generation,reps_formants,reps_params,reps_nishimura) = getValues(elite, n_targets, n_formants, 0, n_var_nishimuras)
                    
                    reps_nishimura = [sum([u,v[-n_global_nishimuras:]],[]) for (u,v) in zip(reps_nishimura,targets_nishimura)] 
                    data = zip(vowels,targets_formants, targets_params, targets_nishimura, alts_formants, alts_params, alts_nishimura, reps_formants, reps_params, reps_nishimura)
                    
                    # wrtite new file
                    for new_line in data:
                        vowel = new_line[0]
                        
                        try:
                            i_rep = int(rep[3:])
                        except ValueError:
                            rep_gen = rep.split(".")
                            (i_rep,i_gen) = (int(rep_gen[0][3:]),int(rep_gen[1]))
                        
                        line_output = [set[:set.index(".")],vowel,i_rep,i_gen]
                        
                        for line_group in new_line[1:]:
                            line_output += [float(v) for v in line_group]
                            
                        line_output.append(int(float(generation)))
                        
                        output.append(line_output)
                    
                    # generate XML
                    root_node = tree.getroot()
                    node_vt = root_node.findall("vocal_tract_model")[0]
                    node_shapes = node_vt.findall("shapes")[0]                   
                    
                    if irep == 0:
                        reps_params = [targets_params[0]] + reps_params
                        target_header = ["target." + v for v in header_params[0]]
                        header_params = [target_header] + header_params
                    
                    for (params_sound_rep,header_sound_rep) in zip(reps_params,header_params):                        
                        vowel = header_sound_rep[0][:header_sound_rep[0].index("_")]
                        name = set[:set.index(".")] + "_" + rep + "_" + vowel
                        header_sound_rep = [item[item.index("_")+1:] for item in header_sound_rep]
                        
                        writeXML(tree, name, header_fixed, params_fixed, header_sound_rep, params_sound_rep)   
    
    print zip(sets,n_maxed_out)
    print sum(n_maxed_out)
    tree.write(os.path.join(root, "JD2.speaker"))
    
    with open(os.path.join(root,"summary.csv"), 'wb') as output_file:
        writer = csv.writer(output_file)    
        for line in output:
            writer.writerow(line)
