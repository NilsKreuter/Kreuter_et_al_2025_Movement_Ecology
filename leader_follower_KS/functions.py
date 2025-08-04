import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy

def ks_2samp(data1, data2):
    """
    Computes the Kolmogorov-Smirnov statistic on 2 samples.
    This is a two-sided test for the null hypothesis that 2 independent samples
    are drawn from the same continuous distribution. It is an asymetric version.
    Parameters
    ----------
    a, b : sequence of 1-D ndarrays
        two arrays of sample observations assumed to be drawn from a continuous
        distribution, sample sizes can be different
    Returns
    -------
    d, D : float
        KS statistic
    prob, p-value : float
        two-tailed p-value
    tau : same type as data
        value of data at which the two cumulative distributions have larger
        difference
    Notes
    -----
    This tests whether 2 samples are drawn from the same distribution. Note
    that, like in the case of the one-sample K-S test, the distribution is
    assumed to be continuous.
    This is the two-sided test, one-sided tests are not implemented.
    The test uses the two-sided asymptotic Kolmogorov-Smirnov distribution.
    If the K-S statistic is small or the p-value is high, then we cannot
    reject the hypothesis that the distributions of the two samples
    are the same.
    Examples
    --------
    >>> from scipy import stats
    >>> np.random.seed(12345678)  #fix random seed to get the same result
    >>> n1 = 200  # size of first sample
    >>> n2 = 300  # size of second sample
    For a different distribution, we can reject the null hypothesis since the
    pvalue is below 1%:
    >>> rvs1 = stats.norm.rvs(size=n1, loc=0., scale=1)
    >>> rvs2 = stats.norm.rvs(size=n2, loc=0.5, scale=1.5)
    >>> stats.ks_2samp(rvs1, rvs2)
    (0.20833333333333337, 4.6674975515806989e-005)
    For a slightly different distribution, we cannot reject the null hypothesis
    at a 10% or lower alpha since the p-value at 0.144 is higher than 10%
    >>> rvs3 = stats.norm.rvs(size=n2, loc=0.01, scale=1.0)
    >>> stats.ks_2samp(rvs1, rvs3)
    (0.10333333333333333, 0.14498781825751686)
    For an identical distribution, we cannot reject the null hypothesis since
    the p-value is high, 41%:
    >>> rvs4 = stats.norm.rvs(size=n2, loc=0.0, scale=1.0)
    >>> stats.ks_2samp(rvs1, rvs4)
    (0.07999999999999996, 0.41126949729859719)
    """
    from scipy import asarray
    from scipy.stats import kstwobign
    data1, data2 = map(asarray, (data1, data2))
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    n1 = len(data1)
    n2 = len(data2)
    data1 = np.sort(data1)
    data2 = np.sort(data2)
    data_all = np.concatenate([data1, data2])
    cdf1 = np.searchsorted(data1, data_all, side = 'right') / (1.0*n1)
    cdf2 = np.searchsorted(data2, data_all, side = 'right') / (1.0*n2)
    tau=0
    darray = cdf1 - cdf2
    d = np.max(np.absolute(darray))
    # Note: d signed distance
    if d == -np.min(darray):
        d = -d
        jamfri = np.min(np.where(darray == np.min(darray))[0])
    else:
        jamfri = np.min(np.where(darray == np.max(darray))[0])
    tau = data_all[jamfri]
    en = np.sqrt(n1*n2/float(n1+n2))
    try:
        prob = kstwobign.sf((en + 0.12 + 0.11 / en) * d)
    except:
        prob = 1.0
    return d, prob, tau

def randomize_times(times, ids = []):
    """
    Randomize the times of the point events of all the ids that are given. This
    just reshuffles the event times among all the individuals taking into
    account.
    Parameters
    ----------
    times : dictionary of lists
        The dictionary contains for each element their times of events in a list
    ids : list of ids
        If not given, the reshuffling is global, if some ids are given,
        only those will be used for the reshuffling.
    Returns
    -------
    times_random : dictionary of lists
        For each element a list of reshuffled event times

    """
    from random import shuffle
    times_random = dict()
    if len(ids) == 0:
        ids = list(times.keys())
    Nevents = dict()
    aux = 0
    tlist = []
    N = len(ids)
    for i in range(N):
        idn = ids[i]
        aux += len(times[idn])
        Nevents[idn] = aux
        tlist.extend(times[idn])
    shuffle(tlist)
    aux=0
    for i in range(N):
        idn = ids[i]
        times_random[idn] = tlist[aux:Nevents[idn]]
        aux = Nevents[idn]
        times_random[idn].sort()
    return times_random

def randomize_ietimes(times, ids = []):
    """
    Randomize the times of the point events of all the ids that are given.
    This randomization keeps the starting time of each individual and reshuffles
    its own interevent times.
    Parameters
    ----------
    times : dictionary of lists
        The dictionary contains for each element their times of events in a list
    ids : list of ids
        If not given, the reshuffling is global, if some ids are given,
        only those will be used for the reshuffling.
    Returns
    -------
    times_random : dictionary of lists
        For each element a list of reshuffled event times

    """
    from random import shuffle
    times_random = dict()
    if len(ids) == 0:
        ids = times.keys()
    for idn in ids:
        Nevents = len(times[idn])
        ietlist = [times[idn][i+1]-times[idn][i] for i in range(Nevents-1)]
        shuffle(ietlist)
        t0 = times[idn][0]
        times_random[idn] = [t0]
        for i in range(Nevents-1):
            t0 += ietlist[i]
            times_random[idn].append(t0)
    return times_random



def waiting_times(times, 
                  ids, 
                  tfloat=True, 
                  t_max=0.0):
    """
    Get the waiting times for two individuals
    Parameters
    ----------
    times : dictionary of lists
        The dictionary contains for each element their times of events in a list
    ids : 2 ids for the reshuffling in a list
        If not given, the reshuffling is global, if some ids are given,
        only those will be used for the reshuffling.
    tfloat : boolean variable
        If True the times are taken as floats, if False event times are datetime
        type
    t_max : float
        If different from 0 then only times 
    Returns
    -------
    tab, tba : lists of time differences

    """
    flag = 0
    tab = list()
    tba = list()
    idi = ids[0]
    idj = ids[1]
    imin = min(times[idi])
    jmin = min(times[idj])
    if jmin > imin:
        a = idj
        b = idi
        flag = 1
    else:
        a = idi
        b = idj
        flag = 0
    Na = len(times[a])
    Nb = len(times[b])
    bx = 0
    ax = 0
    if tfloat:
        while ax < Na-1 and bx < Nb-1:
            while times[b][bx] <= times[a][ax] and bx < Nb-1:
                bx += 1
            if bx != Nb-1:
                aux = times[a][ax] - times[b][bx-1]
                dtab = aux
                tab.append(dtab)
            while times[a][ax] <= times[b][bx] and ax < Na-1:
                ax += 1
            if ax != Na-1:
                aux = times[b][bx] - times[a][ax-1]
                dtba = aux
                tba.append(dtba)
    else:
        while ax < Na-1 and bx < Nb-1:
            while times[b][bx] <= times[a][ax] and bx < Nb-1:
                bx += 1
            if bx != Nb-1:
                aux = times[a][ax] - times[b][bx-1]
                dtab = aux.days*24.0*60.0 + aux.seconds/60.0
                tab.append(dtab)
            while times[a][ax] <= times[b][bx] and ax < Na-1:
                ax += 1
            if ax != Na-1:
                aux = times[b][bx] - times[a][ax-1]
                dtba = aux.days*24.0*60.0 + aux.seconds/60.0
                tba.append(dtba)
    tba = list(filter(lambda x: x != 0.0, tba))
    tab = list(filter(lambda x: x != 0.0, tab))
    if flag == 0:
        return tab, tba
    else:
        return tba, tab


def leadership_network(times,location,species_length_filter,plot_path,
                           scheme = 'global',
                           pmax = 1.0,
                           Nruns = 100,
                           min_int = 50,
                           tfloat = True,
                           rand = 't'
                           ):
    """
    Gives back the network of follower-followees with a maximum p-value pmax,
    following a particular reshuffling scheme.
    Parameters
    ----------
    times : dictionary of lists
        The dictionary contains for each element their times of events in a list
    scheme : string
        'global' for a global reshuffling scheme
        'local' for a local reshuffling scheme
    pmax : float (optional)
        maximum p-value allowed for each edge
    Nruns : integer (optional)
        Number of reshufflings used for getting the p-value
    min_int : integer
        minimum number of interactions (waiting times)
    tfloat : boolean variable
        If True the times are taken as floats, if False event times are datetime
        type
    rand : string
        't' reshuffles the event times among all the individuals
        'iet' reshuffles the interevents for each individual
    Returns
    -------
    g : Networkx DiGraph
        Graph containing the information about the follower-followee network.
        The edges have properties such as A_KS, p and tau.
    """
    if scheme == 'global':
        return A_KS_tau_pvalue_global(times,location,species_length_filter,plot_path,
                                   pmax = pmax,
                                   Nruns = Nruns,
                                   min_int = min_int,
                                   tfloat = tfloat,
                                   rand = rand)
    
    else:
        raise ValueError("Choose a proper reshuffling scheme:'global'")



def A_KS_tau_pvalue_global(times,location,species_length_filter,plot_path,
                           pmax = 1.0,
                           Nruns = 100,
                           min_int = 50,
                           tfloat = True,
                           rand = 't'
                           ):
    """
    Gives back the network of follower-followees with a maximum p-value pmax,
    following a global reshuffling scheme.
    Parameters
    ----------
    times : dictionary of lists
        The dictionary contains for each element their times of events in a list
    pmax : float (optional)
        maximum p-value allowed for each edge
    Nruns : integer (optional)
        Number of reshufflings used for getting the p-value
    min_int : integer
        minimum number of interactions (waiting times)
    tfloat : boolean variable
        If True the times are taken as floats, if False event times are datetime
        type
    rand : string
        't' reshuffles the event times among all the individuals
        'iet' reshuffles the interevents for each individual
    Returns
    -------
    g : Networkx DiGraph
        Graph containing the information about the follower-followee network.
        The edges have properties such as A_KS, p and tau.
    """

# change times with grouped
    import networkx as nx
    import datetime
    g=nx.DiGraph()
    tlist = []
    for key in times.keys():
        tlist.extend(times[key])
    ids = list(times.keys())
    N = len(ids)
    #new dictionary to save the reshuffled AKS value for plotting
    all_AKS = {}
    global all_interactions
    all_interactions = []
    #If this counter goes up it means there has been AKS values under the given premises found. It then enables the reshufflings and randomisation
    #It it stays at 0 it means that no interactions have been found and therefore reshuffling the other data makes no sense
    counter_for_iterations = 0

            # len of ids
    for i in range(N-1):
        for j in range(i+1,N):
            
            #Calcualte the waiting times for A following B and B following A
            if len(times[ids[i]]) >0 and len(times[ids[j]])>0:
                tab, tba = waiting_times(times, [ids[i], ids[j]], tfloat=tfloat)
            
                if (len(tab) + len(tba)) > min_int:
                    all_interactions += [(len(tab) + len(tba))]
            else:
                tab = []
                tba = []
            
	    #AKS values are only calculated if the number of interactions previously calculated are more than the ones required. 
            if len(tab)> 0 and len (tba) > 0 and (len(tab) + len(tba)) > min_int:

                A_KS, p_bad, tau = ks_2samp(tab, tba)
                
                interaction_counter = len(tab) + len(tba)
               	
                
                #Original timestamps
                ind_one = pd.DataFrame({'floats': times[ids[i]] })
                ind_one['datetime'] = pd.to_datetime(ind_one['floats'], unit='s')
               
                ind_two = pd.DataFrame({'floats':  times[ids[j]]})
                ind_two['datetime'] = pd.to_datetime(ind_two['floats'], unit='s')

                #Create an empty list for y-values for original data
                y_values_1_org = [1] * len(ind_one)
                y_values_2_org = [2] * len(ind_two)
                
                
                # Plot the 'Timestamps' column with ticks and no y-axis
                plt.figure(figsize=(10, 2), dpi=300)  # Adjust the figure size as needed
                plt.ylim(0.5, 2.5)
                plt.yticks([1,2], [ids[i],ids[j]])
                plt.plot(ind_one['datetime'], y_values_1_org, marker='|', markersize=20, linestyle='None', color='b')
                plt.plot(ind_two['datetime'], y_values_2_org, marker='|', markersize=20, linestyle='None', color='r')
                
                plt.xlabel('Timestamps')
                plt.title(f'Timestamp REAL   interactions: {interaction_counter}  {location} {species_length_filter}')
                plt.grid(axis='x')
                #plt.gca().yaxis.set_major_locator(plt.NullLocator())  # Hide y-axis labels
                
                filename_time = f"{plot_path}{location}_{ids[i]}_{ids[j]}__TIME_plot.png"

                
                plt.savefig(filename_time)
                plt.close()
                #If this counter goes up it means there has been AKS values under the given premises found. It then enables the reshufflings and randomisation
                counter_for_iterations = counter_for_iterations +1
            
            else:
                A_KS, p_bad, tau=(0.0, 0.0, 0.0)
            #Negative AKS are flipped into positive to check again
            
            ##### CHECK this whether it changes the ids as well and is plotted correctly and given out in the correct order!
            
            if A_KS != 0:
                if A_KS < 0.0:
                    g.add_edge(ids[j], ids[i], A_KS = -A_KS, tau=tau, p=Nruns)
                else:
                    g.add_edge(ids[i], ids[j], A_KS = A_KS, tau=tau, p=Nruns)
                

    if counter_for_iterations > 0:
        for irun in range(Nruns):
            
            #print(Nruns-irun)
            #Wahtever randomisation scheme was chosen the times will the randomised
            
            #For each edge (interactions of two individuals) the "true" AKS value of the sample data will be saved for later plotting
            for edge in g.edges():
                i = edge[0]
                j = edge[1]
                
               
                #This is where the AKS value is saved for the "true" value of the sample
                A_KS = g[i][j]['A_KS']
                

                #
                
                #Randomize the data only for the two individuals
                if rand == 't':
                    t_rand = randomize_times(times)
                elif rand == 'iet':
                    t_rand = randomize_ietimes(times)
                

                #then do the waiting times 
                tab, tba = waiting_times(t_rand, [i, j], tfloat = tfloat)
                #The reshuffled AKS value is saved for later plotting. THE AKS rand value was set to 1 in case it does not meet the criterium
                if len(tab)> 0 and len (tba) > 0 and (len(tab) + len(tba)) > min_int:
                    A_KS_rand, p_bad, tau = ks_2samp(tab, tba)
                else:
                    A_KS_rand, p_bad, tau=(1, 0.0, 0.0)
                    
                #The p value was saved as the number of runs before and now it checks if the randomised AKS value is bigger or smaller than the "true" AKS value
                #The AKS value is taken as the absolute value and therefore always positive. Therefore it is possible against the "true" AKS value that will
                #always be positive. 
                #the "p" vlaue that conssist of the Nruns number will be deducted by 1 in case the randomised AKS value is lower than the "true" AKS value
                if abs(A_KS_rand) < abs(A_KS):
                    g[i][j]['p'] -= 1
    
                #save the AKS values for plotting. If the edge is new it will create a new entry in the reshufffled AKS value dicitionary for the first value
                if edge not in all_AKS:
                    all_AKS[edge] =[]
                all_AKS[edge].append(A_KS_rand)
        
                

    
    for edge in g.edges():
        i = edge[0]
        j = edge[1]
        #the p value is being calculate here. the current p" value is divided by the number of runs. This gives the indication of whether the "true"
        #AKS value falls under the distribution of the reshuffled ones
        g[i][j]['p'] = float(g[i][j]['p'])/float(Nruns)
    G=nx.DiGraph( [ (u,v,d) for u,v,d in g.edges(data=True) if d['p'] <= pmax] )
           
    
    #Create new dicitonary for all pairs that had a AKS above 0
    filtered_edges_dict_AKS = {}
    #Create new dicitonary for all pairs that had a P value below 0.05
    filtered_edges_dict_p = {}
    
    for u, v, data in g.edges(data=True):
        AKS_value = data['A_KS']
        p_value = data['p']
        #and p_value < 0.05
        if AKS_value > 0 :
            edge_name = (u,v)
            filtered_edges_dict_AKS[edge_name] = float(AKS_value)
            filtered_edges_dict_p[edge_name] = float(p_value)
    
    
    counter = 0
    #Plot each interaction with AKS above 0
    while counter < len(filtered_edges_dict_AKS.keys()):
        
        sampled_AKS = filtered_edges_dict_AKS[list(filtered_edges_dict_AKS.keys())[counter]]
        sampled_p = filtered_edges_dict_p[list(filtered_edges_dict_p.keys())[counter]]
            
        shuffled_AKS = all_AKS[list(filtered_edges_dict_AKS.keys())[counter]]
    
     
        plt.figure(dpi=300)
        bin_range = (-1, 1)  # Your desired range
        num_bins = 200
        
        hist, bin_edges = np.histogram(shuffled_AKS, bins=num_bins, range=bin_range)
        
        bin_width = (bin_range[1] - bin_range[0]) / num_bins
        bin_centers = bin_edges[:-1] + bin_width / 2
        
        # Plot the histogram using matplotlib with 'step' style
        plt.step(bin_centers, hist, where='mid', color='black', label='Reshuffled')

        plt.xlabel('AKS')
        plt.ylabel('Frequency')
 	#add {location} {species_length_filter}
        plt.title(f'{list(filtered_edges_dict_AKS.keys())[counter]} {location} {species_length_filter}')
        plt.axvline(sampled_AKS, color='red', linestyle='--', label='Real')
        plt.text(sampled_AKS + 0.1, (max(hist)/2), f'p={sampled_p}', rotation=90)
        plt.legend(loc='upper left')
        
        filename_AKS = f"{plot_path}{location}_{list(filtered_edges_dict_AKS.keys())[counter]}_AKS_plot.png"

        
        plt.savefig(filename_AKS)      
        plt.close()
        counter = counter + 1

   
    return G

