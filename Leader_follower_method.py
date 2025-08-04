######################################################################################################################
##################################################### User input  ###################################################
######################################################################################################################

################ Insert the information needed for the script to run #############

#Set the number of CPU cores used for analysis. The script uses multiprocessing to reduce analysis time. (Pay attention to specific requirements from your IDE and available cores.)
cpus_used = 1

#Set the path to your .csv data file (Exp: "C:/Users/Example/data.csv")
csv_path_raw_data = ""

#Set the timezone in which the data is saved
timezone_raw_data = "UTC"

#Set the timezone in which the data was collected
timezone_location = "UTC"

#Specify the species (this variable will be used for the plotting of graph titles and other output files (see GitHub Readme))
species_length_filter = "Tiger"

#Specifically for IMOS ATF data. Mention the length units mentioned in the data. Either "cm", "mm" or "cm and mm". The script will take care of the letters in the length measurements of the IMOS ATF data. 
#The code needs to be updated in case you use different data (Transforming the measurements section)
species_length_unit = "mm"

#Set the maximum p value for which a leadership pattern is determined (see manuscript for more details)
pmax = 0.05

#Set the number of reshufflings for all null-model datasets in every analysis (see manuscript for more details)
Nruns = 10000

#Set the number of interactions that a pair must have to be considered for analysis (see manuscript for more details)
min_int = 200

#Set the output path where the graphs and other output data should be saved (Exp: "C:/Users/Example/output/") . Subfolders and all output will be saved here
output_path = ""


######################################################################################################################
########################################################### Start of script ##########################################
######################################################################################################################

###########################
# Import all packages used 
###########################

import multiprocessing
import leader_follower_KS.functions
import datetime
import pandas as pd
import itertools #
import matplotlib.pyplot as plt
import seaborn as sns
import os 
import random 
import networkx as nx
import re
from scipy import stats
from scipy.stats import linregress
from itertools import combinations
import pytz
import copy
import numpy as np
from networkx import gnm_random_graph, connected_components, descendants, directed_configuration_model
from scipy.stats import gaussian_kde


#########################################################
#Function to create character IDs for numerical tag IDs
#########################################################

def generate_character_ids(num_ids):
    # Generate a sequence of uppercase letters
    letters = itertools.cycle('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

    # Generate one-character IDs until the letters run out
    for _ in range(26):
        yield next(letters)

    # Generate two-character IDs once the letters are depleted
    for i in range(26, num_ids):
        for first_letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for second_letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                yield f"{first_letter}{second_letter}"
                if i == num_ids - 1:
                    return
                i += 1



#############################
# Upload and prepare raw data 
#############################

#Load the raw data
raw_data = pd.read_csv(csv_path_raw_data)

#Change column names for id, timestmap and sex. station names provided by IMOS is not used and renamed. 
columns_name={'tag_id': 'id', 'detection_datetime': 'timestamp', 'animal_sex': 'sex', 'station_name': 'IMOS_station_name'}
raw_data = raw_data.rename(columns=columns_name)

#set the timestamp to the timezone of the raw data. This helps to transform it to the desired timezone of the location data
if (timezone_raw_data == "UTC"): 
    raw_data['timestamp'] = pd.to_datetime(raw_data['timestamp'], utc=True)
else:
    raw_data['timestamp'] = raw_data['timestamp'].dt.tz_localize(pytz.timezone(timezone_raw_data))

#convert the timezone to location timezone
if timezone_location == "UTC" :
    print("Timezone raw data matches location timezone")
else:
    raw_data['timestamp'] = raw_data['timestamp'].dt.tz_convert(pytz.timezone(timezone_location))

raw_data = raw_data.sort_values(by='timestamp')

#################
# Receiver naming 
#################

#There are instances that IMOS stations have the same name and different positions.

#New receiver naming scheme for a new station_name column
#Get all receiver positions based on the longitude and latitude combinations in the data
receiver_information = raw_data.drop_duplicates(subset=['receiver_deployment_longitude', 'receiver_deployment_latitude'])[['receiver_deployment_longitude','receiver_deployment_latitude']]

#Create new names for the speciic receiver positions. They are now called Receiver_1, Receiver_2 etc. to avoid conflicts for same IMOS receiver names at different positions
unique_numbers = {}
receiver_information['station_name'] = receiver_information.apply(lambda row: f"Receiver_{unique_numbers.setdefault((row['receiver_deployment_longitude'], row['receiver_deployment_latitude']), len(unique_numbers) + 1)}", axis=1)

# The combination of longitude and latitude are checked and the new station name is added corresponding to the receiver_positions names
raw_data = pd.merge(raw_data, receiver_information, on=['receiver_deployment_longitude', 'receiver_deployment_latitude'], how='left')


###################################################################################
# Create an overview dataframe for all indiviauls and setting numeric character IDs 
###################################################################################
#Using single or double characters isntead of numeric strings improves indetification of individuals in plots, output etc.

#Get the unique ids and create the unique character IDs for them. This list will later be used to get the correct character ID for the numeric ID
allIds = pd.DataFrame(raw_data['id'].unique(), columns=['ID_numeric'])
allIds['ID_character'] = list(itertools.islice(generate_character_ids(len(allIds)), len(allIds)))

# Get the measuremtens into the allIds, already transformed into total length (cm)
# Create new columns for the ID and the measurement for this ID
allIds = pd.merge(allIds, raw_data[['id','measurement']], left_on='ID_numeric', right_on='id', how='left')
#Drop all duplicates for each ID
allIds.drop_duplicates(subset='ID_numeric', inplace=True)
#Drop the second ID column that was created in this process
allIds.drop(columns=['id'], inplace=True)
#Reset index to ensure better use later
allIds= allIds.reset_index(drop=True)
 

################################
# Transforming the measurements
###############################

#Fork Length to Length conversion. 
#Paramaters for conversion are based on the species input by the user
if species_length_filter =="Grey":
    value1 =1.130
    value2 =7.67 

#Grey Total length in cm: 1.130*Fork Length(cm) + 7.67 
#Bradley, D., Conklin, E., Papastamatiou, Y. P., McCauley, D. J., Pollock, K., Kendall, B. E., ... & Caselle, J. E. (2017). Growth and life history variability of the grey reef shark (Carcharhinus amblyrhynchos) across its range. PLoS One, 12(2), e0172370. Table 1

if species_length_filter =="Blacktip":
    value1 =1.16
    value2 =4.16
    
# Blacktip: Total length in cm : 1.16*Fork Length(cm) + 4.16 
#Lyle, J. M. (1987). Observations on the biology of Carcharhinus cautus (Whitley), C. melanopterus (Quoy & Gaimard) and C. fitzroyensis (Whitley) from northern Australia. Marine and Freshwater Research, 38(6), 701-710.
    
if species_length_filter =="Tiger":
    value1 =1.12
    value2 =11.74

# Tiger: Total length in cm: 1.12*Fork Length(cm) + 11.74 
#Bartes, S. N., & Braccini, M. (2023). Length length relationships for the main shark species caught in the commercial shark fisheries of Western Australia. Fisheries Management and Ecology, 30(2), 224-227.

#Define the patter that is used for the transformation of the measurements. Used to remove the characters and find the correct numbers in the measurement column.
pattern = r"\b\d+\b"

#################################################################################################################
# Transforming the measurements into cm and total Length. Fork Length is calculated based on the defintions above
#################################################################################################################

if species_length_unit == "cm and mm":
    
    #extract the Fork length values and the non Fork lengths. Non Fork Lengths are all treated as Total Length
    no_length_calculation = allIds[~allIds['measurement'].str.contains('FORK LENGTH', case=False,na=True)].reset_index(drop=True)
    length_calculation = allIds[allIds['measurement'].str.contains('FORK LENGTH', case=False, na=False)].reset_index(drop=True)

    #Extract the Ids that need a transformation from fork length to total length in both mm and cm
    length_calculation_mm = length_calculation[length_calculation['measurement'].str.contains('mm', case=False,na=False)].reset_index(drop=True)
    length_calculation_cm = length_calculation[length_calculation['measurement'].str.contains('cm', case=False,na=False)].reset_index(drop=True)


########### MM #############   
    #Check if measuremnts in mm were found
    if len(length_calculation_mm)>0:
        
        #If species length unit is mm the values are calcualted based on the values extracted from the publications and divided by 10 to get cm unit values
        length_calculation_mm['measurement_new_total_cm'] = length_calculation_mm['measurement'].apply(lambda x: (int(re.search(pattern, x).group())/10 * value1 + value2) if re.search(pattern, x) else None)
        
########### CM #############   
    #Check if measuremnts in mm were found
    if len(length_calculation_cm)>0:

        #If species length unit is mm the values are calcualted based on the values extracted from the publications and divided by 10 to get cm unit values
        length_calculation_mm['measurement_new_total_cm'] = length_calculation_mm['measurement'].apply(lambda x: (int(re.search(pattern, x).group()) * value1 + value2) if re.search(pattern, x) else None)
        

    #Fuse both dataframes that were calculated from Fork Length to Total length together again
    length_calculation_complete =pd.concat([length_calculation_mm, length_calculation_cm], ignore_index=True)
    
################### remove the characters for the ones in total length and transform some values from mm to cm if needed ###### 
    #Extract the Ids that need a transformation from mm in cm and remove the characters for the ones with cm
    correct_length_mm = no_length_calculation[no_length_calculation['measurement'].str.contains('mm', case=False,na=False)].reset_index(drop=True)
    correct_length_cm = no_length_calculation[no_length_calculation['measurement'].str.contains('cm', case=False,na=False)].reset_index(drop=True)
    
########### MM #############   
    correct_length_mm['measurement_new_total_cm'] = correct_length_mm['measurement'].apply(lambda x: (int(re.search(pattern, x).group())/10)) 

########### CM #############   
    correct_length_cm['measurement_new_total_cm'] = correct_length_cm['measurement'].apply(lambda x: (int(re.search(pattern, x).group()))) 

#Fuse both dataframes that were calculated from Fork Length to Total length together again
    all_calculations =pd.concat([correct_length_mm, correct_length_cm,length_calculation_complete], ignore_index=True)
    all_calculations['measurement_new_total_cm'] = all_calculations['measurement_new_total_cm'].astype(int)

###################################################################################
# Conversion of species length into total length inf length unit is either mm or cm
###################################################################################

if species_length_unit == "mm" or species_length_unit == "cm":

    #Two dataframes are created to separate the data into length measurements that need to be converted (length_calculation) and the ones that don't need to be converted (no_length_calcualtion) but need to have their characters removed in the column
    #extract the Fork length values and the non Fork lengths. Non Fork Lengths are all treated as Total Length
    no_length_calculation = allIds[~allIds['measurement'].str.contains('FORK LENGTH', case=False,na=True)].reset_index(drop=True)
    length_calculation = allIds[allIds['measurement'].str.contains('FORK LENGTH', case=False,na=False)].reset_index(drop=True)
    
  
    ########### MM #############

    #If species length unit is mm the values are calcualted based on the values extracted from the publications and divided by 10 to get cm unit values
    if species_length_unit == "mm":        
        length_calculation['measurement_new_total_cm'] = length_calculation['measurement'].apply(lambda x: (int(re.search(pattern, x).group())/10 * value1 + value2) if re.search(pattern, x) else None)
    
    #For the individuals that are already in total length this part removes the characters in the column and divides the value by 10 to get cm unit values
    if species_length_unit == "mm":            
        no_length_calculation['measurement_new_total_cm'] = no_length_calculation['measurement'].apply(lambda x: (int(re.search(pattern, x).group())/10)) 

    ########### CM #############

    #If species length unit is cm the values are calcualted based on the values extracted from the publications
    if species_length_unit == "cm":
        length_calculation['measurement_new_total_cm'] = length_calculation['measurement'].apply(lambda x: (int(re.search(pattern, x).group()) * value1 + value2) if re.search(pattern, x) else None)
    
    #For the individuals that are already in total length this part removes the characters in the column to get cm unit values
    if species_length_unit == "cm":
        no_length_calculation['measurement_new_total_cm'] = no_length_calculation['measurement'].apply(lambda x: (int(re.search(pattern, x).group()))) 

    #Fuse the no_length and length dataframes together 
    all_calculations =pd.concat([no_length_calculation, length_calculation], ignore_index=True)
    all_calculations['measurement_new_total_cm'] = all_calculations['measurement_new_total_cm'].astype(int)


################################################################################
# Fuse the new calculations into the allIds object and omit the old measurements
################################################################################

#Remove the old measurement column from the allids dataframe
allIds.drop(columns='measurement', inplace=True)

#Merge the newly calculated measurements with the allids object    
allIds = allIds.merge(all_calculations[['ID_numeric', 'measurement_new_total_cm']], on='ID_numeric', how='left')

#rename the measuremnt column to indicate that the cm unit
allIds.rename(columns={'measurement_new_total_cm': 'measurement_cm'}, inplace=True)

###############################
# Insert character ID into raw 
###############################

#rename the raw data measuremnt column that uses differnt units and lenght measuremnts. It will not be used anymore
raw_data.rename(columns={'measurement': 'measurement_old'}, inplace=True)

raw_data = pd.merge(raw_data, allIds, left_on='id', right_on='ID_numeric')

#########################
## Insert sex into allIds 
#########################

allIds = allIds.merge(raw_data[['ID_numeric', 'sex']].drop_duplicates(subset=['ID_numeric']), on='ID_numeric', how='left')

########################################
# Prepare data for multiprocessing steps 
########################################
#All data will be saved in a nested dict that the multiprocessing can use

#Save all receiver names
receiver_names = raw_data['station_name'].unique()

#set a counter that is used to cycle through the receiver names in the following loop
receiver_count = 0

#create an empty list and empty dataframe to get the information for leader and follower for each individual at each receiver

LF_info_final = pd.DataFrame()

all_data_dict = {}


while receiver_count < len(receiver_names):

      
    # Select a specific location of receivers (especially when using the leader follower network)
    location = receiver_names[receiver_count]
        
    #Filter raw data for set location
    times_location = raw_data[raw_data['station_name'] == location]
    
    #Use that filtered data and extract IDs and timestamps 
    times_ID_TIMESTAMP = times_location[['ID_character','timestamp']]

    #transform data into floats
    timestamps = times_ID_TIMESTAMP['timestamp']
    timestamp_floats = []
    timestamp_floats = timestamps.view('int64') / 10**9
    
    #Create separate object with ID
    ind_all = times_ID_TIMESTAMP['ID_character']

    #Merge ID and JD timestamps
    times_ID_TIMESTAMP =  pd.DataFrame({"ID_character": ind_all, "timestamp": timestamp_floats})
    
    #transforming the data for the next step
    grouped = times_ID_TIMESTAMP.groupby('ID_character')['timestamp'].apply(list)

    #convert the resulting Series to a dictionary
    grouped = grouped.to_dict()

    all_data_dict[location] = grouped

    receiver_count = receiver_count+1



##########################################
#Organising data for Multicore Processing
##########################################

#Convert your main dictionary (all_data_dict) into a list of tuples. Each tuple will contain the receiver's name and the corresponding dictionary. This structure makes it easier to distribute the workload across multiple cores.
data_list = [(name, all_data_dict[name]) for name in all_data_dict]


#########################################
## Worker function for multicore process 
########################################
# Worker Function: This function will contain the code that processes each dictionary. It should accept a tuple as an argument and return the result of your analysis.

def process_dictionary(data_dict_all,plot_path):
    processed_data = []
    try:
        #Due to the multicore processing the function will take each receiver data separately
        receiver_name, data_dict= data_dict_all
        print(receiver_name)
        
        #main function that is called from the leader_follower_KS package
        g = leader_follower_KS.functions.leadership_network(data_dict,receiver_name,species_length_filter,plot_path,
                                                       scheme = 'global',
                                                       pmax = pmax,
                                                       Nruns = Nruns,
                                                       min_int = min_int,
                                                       tfloat = True,
                                                       rand = 'iet')
        
        #If significant interactions are found among individauls at the specific receiver the plotting the graphs
        
        if g.number_of_nodes() > 0:
            #Filter all the Ids that are represented in the nodes
            allIds_filter = allIds[allIds['ID_character'].isin(list(g.nodes))]

            #Add the sex as metadata to the network
            for index, row in allIds_filter.iterrows():
                node_id = row['ID_character']
                node_attributes = {'sex': row['sex']}
                g.add_node(node_id, **node_attributes)

            #Add the size into the metadata. Size in cm will be used as node size
            for index, row in allIds_filter.iterrows():
                node_id = row['ID_character']
                node_attributes = {'measurement': (row['measurement_cm'])}
                g.add_node(node_id, **node_attributes)

            #create node sizes for the network
            node_sizes = 1000

            ### get the number of incoming and outgoing nodes
            # Initialize dictionaries to store incoming and outgoing degrees
            incoming_degrees = {}
            outgoing_degrees = {}

            # Iterate through nodes
            for node in g.nodes():
              incoming_degrees[node] = sum(1 for _ in g.predecessors(node))
              outgoing_degrees[node] = sum(1 for _ in g.successors(node))

            #Create the node labels that are plotted. This info is also used to be stored for each recevier to check the raw numbers at the end
            node_labels = {node: f"{node}\n L:{incoming_degrees} F:{outgoing_degrees}" for node, (incoming_degrees, outgoing_degrees) in
                      zip(g.nodes, zip(incoming_degrees.values(), outgoing_degrees.values()))}
            
            #This is splitting up the node_labels into a dataframe
            for key, value in node_labels.items():
                parts = value.split('\n')
                name = parts[0]
                values = parts[1].split()
                l_value = int(values[0].split(':')[1])
                f_value = int(values[1].split(':')[1])
                processed_data.append({'Name': name, 'L': l_value, 'F': f_value, 'location': receiver_name})

            # Creating a DataFrame from the processed data
            LF_info_temp = pd.DataFrame(processed_data)
            
          
            # Initialise an empty list to store edge strength of the data
            LF_strength_info_list = []
            
            # Iterate over all edges in the graph
            for edge in g.edges(data=True):
                # edge is a tuple in the form of (node1, node2, data_dict)
                node1, node2, data = edge
                # Extract the 'AKS' value from the edge's attribute dictionary
                AKS_value = data.get('A_KS', None)  

                # Get the 'sex' attribute for each node
                sex_node1 = g.nodes[node1].get('sex', None)
                sex_node2 = g.nodes[node2].get('sex', None)
                
                # Get the 'measurement' attribute for each node
                measurement_node1 = g.nodes[node1].get('measurement', None)
                measurement_node2 = g.nodes[node2].get('measurement', None)
                
                # Append an entry for the source node (node1) with behavior "Follower" and target node as "cpd_ID"
                LF_strength_info_list.append((node1,sex_node1,measurement_node1, "Following", -AKS_value, node2))
                # Append an entry for the target node (node2) with behavior "Leader" and source node as "cpd_ID"
                LF_strength_info_list.append((node2,sex_node2,measurement_node2, "Leading", AKS_value, node1))
                
            # Create a DataFrame from the LF_strength_info_temp list
            LF_strength_info_temp = pd.DataFrame(LF_strength_info_list, columns=['ID', 'Sex', 'Length', 'Behaviour', 'AKS', 'cpd_ID'])
            LF_strength_info_temp['AKS'] = pd.to_numeric(LF_strength_info_temp['AKS'])
            LF_strength_info_temp['Location'] = receiver_name
            
            #################
            #plot the network
            #################
            pos = nx.shell_layout(g)  # You can choose different layout algorithms
            
            # Separate nodes based on sex
            male_nodes = [node for node in g.nodes if g.nodes[node]['sex'] == 'MALE']
            female_nodes = [node for node in g.nodes if g.nodes[node]['sex'] == 'FEMALE']
                  
            # Calculate padding for the plot limits
            x_values, y_values = zip(*pos.values())
            x_min, x_max = min(x_values), max(x_values)
            y_min, y_max = min(y_values), max(y_values)
            padding = 0.35  # Adjust this value to change the amount of padding
            
            plt.figure(dpi=300)
            plt.title(f'{species_length_filter} - {receiver_name}')
            
            # Set plot limits with padding
            plt.xlim(x_min - padding, x_max + padding)
            plt.ylim(y_min - padding, y_max + padding)
            
            # # Draw edges and increasing the node size to reduce overlapping from nodes and edges
            nx.draw_networkx_edges(g, pos, edge_color="grey", arrows=True, arrowsize=10, node_size=node_sizes*2, 
                                   width=[data.get('A_KS', 1)*35 for _, _, data in g.edges(data=True)])
            
            # Draw male nodes as squares with white fill and black edges
            nx.draw_networkx_nodes(
                g, pos, nodelist=male_nodes, node_shape='s', node_color='darkgrey',
                edgecolors='black', node_size=node_sizes
            )
            
            # Draw female nodes as triangles with white fill and black edges
            nx.draw_networkx_nodes(
                g, pos, nodelist=female_nodes, node_shape='D', node_color='white',
                edgecolors='black', node_size=node_sizes
            )
                
                       
            #Labelling the network nodes
            nx.draw_networkx_labels(g, pos,labels={key: key for key in node_labels}, font_size=18,font_color="black",  verticalalignment='center',
                                   horizontalalignment='center')


             # Create the legend for shapes
            legend_elements = [
                plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='darkgrey', 
                           markeredgecolor='black', markersize=10, label='Male'),
                plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='white', 
                           markeredgecolor='black', markersize=10, label='Female')
            ]
            
            plt.legend(handles=legend_elements, bbox_to_anchor=(1, 1),frameon=True, edgecolor='black', loc='upper left')
            # plt.show()
            
            #######################################################
            ##creating the path and filename for saving the network
            #######################################################
            filename_network = f"{plot_path}Network_plot_{receiver_name}_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
            plt.tight_layout()
            plt.savefig(filename_network)


            plt.close()
            
            
        else:
            LF_info_temp = pd.DataFrame() 
            LF_strength_info_temp = pd.DataFrame()



    
        return LF_info_temp,LF_strength_info_temp,g, {
    'receiver': receiver_name,
    'edges': g.number_of_edges(),
    'nodes': g.number_of_nodes(),
    'species': species_length_filter,
    'min_interaction': min_int
}
    except Exception as e:
        print(f"Error processing {data_dict_all[0]}: {e}") 
        return pd.DataFrame() 
    
####################################
#Distributing the Work Across Cores 
####################################

if __name__ == "__main__":
    
    ##################################################
    # Create output paths for the graphs and csv files
    ##################################################
    timestamp = datetime.datetime.now()
    time_date_output = timestamp.strftime("%Y%m%d_%H%M%S")
    print(time_date_output)
    plot_path = f'{output_path}{species_length_filter}_{time_date_output}/plots/'
    csv_path = f'{output_path}{species_length_filter}_{time_date_output}/'

    # Create the directories using os.makedirs
    os.makedirs(plot_path, exist_ok=True)
    os.makedirs(csv_path, exist_ok=True)
    
    
    #########################################
    # Plot all detections for all individuals
    #########################################
    
    # Group by ID_character and find the first and last detection
    detection_times = raw_data.groupby('ID_character')['timestamp'].agg(['min', 'max']).reset_index()

    # Plot the timeline
    fig, ax = plt.subplots(figsize=(10, 5), dpi =300)

    for idx, row in detection_times.iterrows():
        ax.plot([row['min'], row['max']], [row['ID_character'], row['ID_character']], marker='o')

    ax.set_ylabel('ID')
    ax.set_title(f'First and last detection ({species_length_filter})')
    plt.xticks(rotation=45)
    plt.tight_layout()


    filename_detection_plot = f"{plot_path}detection_plot_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    plt.savefig(filename_detection_plot)
    plt.close()


    with multiprocessing.Pool(processes=cpus_used) as pool:
        print("Starting multiprocessing pool")
        # Collect all returned dataframes in a list
        results = pool.starmap(process_dictionary,  [(data, plot_path) for data in data_list])
        print("Multiprocessing pool finished")
        
        
    ########################################
    # Analysing the output for all receivers 
    ########################################
    
     ## Unpack results into DataFrames and graphs
    list_of_dfs,list_of_strength, list_of_graphs, additional_info  = zip(*results)
    
   #Save the interaction test output
    additional_info_df = pd.DataFrame(list(additional_info))
    additional_info_df.to_csv(f'{csv_path}interaction_output_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.csv', index=False)
   
     # Concatenate all dataframes in the list
    final_LF_info = pd.concat(list_of_dfs, ignore_index=True)
    
    # Concatenate all dataframes in the list
    final_list_of_strength = pd.concat(list_of_strength, ignore_index=True)
    final_list_of_strength.to_csv(f'{csv_path}LF_strength_{species_length_filter}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.csv',index=False)
    
    #save information of all IDs
    allIds.to_csv(f'{csv_path}all_ID_{species_length_filter}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.csv', index=False)

    #save the information for all receivers
    receiver_information.to_csv(f'{csv_path}receiver_information_{species_length_filter}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.csv', index=False)
    
    
    ####################################
    # Plotting grouped receiver networks 
    ####################################
        
    if len(final_list_of_strength)>0: 
        
        
        ########## Prepare data for scatterplot and boxplot ##############
        summed_values = final_list_of_strength.groupby('ID')['AKS'].sum().reset_index()

        summed_values = final_list_of_strength.groupby('ID').agg({
            'AKS': 'sum',       # Sum the AKS values
            'Sex': 'first',     # Keep the first non-null sex value
            'Length': 'first'   # Keep the first non-null length value
        }).reset_index()


        # Separate the data by sex
        female_data = summed_values[summed_values['Sex'] == 'FEMALE']
        male_data = summed_values[summed_values['Sex'] == 'MALE']
        
        #exclude IDs for scatterplot that don't have length measurements in the data. sex specific and all ids together for the overall trend
        summed_values_length_recorded = summed_values.dropna(subset=['Length'])
        female_data_length_recorded = female_data.dropna(subset=['Length'])
        male_data_length_recorded = male_data.dropna(subset=['Length'])
        
        #exclude IDs for boxplots that don't have sex measurements in the data
        summed_values_sex_recorded = summed_values.dropna(subset=['Sex'])


        ############################
        # SCATTER and BOXPLOT OUTPUT 
        ############################
        
        #Create scatterplot when bot sexes are found and only create the boxplot if both are found

        if len(female_data_length_recorded)>0 and len(male_data_length_recorded)>0:
            
            # Perform linear regression for combined data
            slope_all, intercept_all, r_value_all, p_value_all, std_err_all = linregress(summed_values_length_recorded['AKS'], summed_values_length_recorded['Length'])
            
            # Perform linear regression for female data
            slope_female, intercept_female, r_value_female, p_value_female, std_err_female = linregress(female_data_length_recorded['AKS'], female_data_length_recorded['Length'])
            
            # Perform linear regression for male data
            slope_male, intercept_male, r_value_male, p_value_male, std_err_male = linregress(male_data_length_recorded['AKS'], male_data_length_recorded['Length'])
            
            # Plot scatter plot for female and male data
            plt.figure(dpi=300)
            plt.scatter(female_data_length_recorded['AKS'], female_data_length_recorded['Length'], edgecolor='black', facecolor='white', marker='D' , label='Female')
            plt.scatter(male_data_length_recorded['AKS'], male_data_length_recorded['Length'], edgecolor='black', facecolor='darkgrey', marker='s', label='Male')
            
            # Plot trend line for combined data
            plt.plot(summed_values_length_recorded['AKS'], slope_all*summed_values_length_recorded['AKS'] + intercept_all, color='black',linestyle='solid', label='Combined')
            
            # Plot trend line for female data
            plt.plot(female_data_length_recorded['AKS'], slope_female*female_data_length_recorded['AKS'] + intercept_female, color='lightgrey',linestyle='solid', label='Female')
            
            # Plot trend line for male data
            plt.plot(male_data_length_recorded['AKS'], slope_male*male_data_length_recorded['AKS'] + intercept_male, color='darkgrey',linestyle='solid', label='Male')
            
            # Annotate the plot with the equation of the trend line and p-value for combined data
            plt.annotate(f'Combined: y = {slope_all:.2f}x + {intercept_all:.2f} p = {p_value_all:.4f}', 
                         xy=(0.5, (plt.gca().xaxis.label.get_position()[1]-0.22)), xycoords='axes fraction',
                         fontsize=10, ha='center', va='top',
                         bbox=dict(boxstyle='round', fc='white', alpha=0.5))
            
            # Annotate the plot with the equation of the trend line and p-value for female data
            plt.annotate(f'Female: y = {slope_female:.2f}x + {intercept_female:.2f} p = {p_value_female:.4f}', 
                         xy=(0.5, (plt.gca().xaxis.label.get_position()[1]-0.32)), xycoords='axes fraction',
                         fontsize=10, ha='center', va='top',
                         bbox=dict(boxstyle='round', fc='white', alpha=0.5))
            
            # Annotate the plot with the equation of the trend line and p-value for male data
            plt.annotate(f'Male: y = {slope_male:.2f}x + {intercept_male:.2f} p = {p_value_male:.4f}', 
                         xy=(0.5, (plt.gca().xaxis.label.get_position()[1]-0.42)), xycoords='axes fraction',
                         fontsize=10, ha='center', va='top',
                         bbox=dict(boxstyle='round', fc='white', alpha=0.5))
            
            plt.title('Influence of size on leader-follower preference')
            
            plt.xlabel('AKS')
            plt.ylabel('Size(cm)')
            # Place the legend outside of the plot on the right side
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1),frameon=True, edgecolor='black',)
            plt.tight_layout() 
            filename_scatter = f"{plot_path}scatter_plot_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
            plt.savefig(filename_scatter)
            plt.close()
            
            
            ######################## BOXPLOT OUTPUT for Sex and AKS  ################################          
            plt.figure(dpi=300)
            sns.boxplot(x='Sex', y='AKS', data=summed_values_sex_recorded, palette={"FEMALE": "white", "MALE": "darkgrey"})
            
            # Set the y-axis limit to be higher than the max to be able to include a signficiance or non-significance statement
            plt.ylim(summed_values_sex_recorded['AKS'].min()-0.05, summed_values_sex_recorded['AKS'].max()+0.15)  
            
            # Add title and labels
            plt.title('Influence of sex on leader-follower preference')
            plt.xlabel('Sex')
            plt.ylabel('AKS')
            
            # Perform statistical test (Mann-Whitney U test)
            female_AKS = summed_values_sex_recorded[summed_values_sex_recorded['Sex'] == 'FEMALE']['AKS']
            male_AKS = summed_values_sex_recorded[summed_values_sex_recorded['Sex'] == 'MALE']['AKS']
            
            ############ test for the correct stats test ########### 
            # Normality tests
            normality_female = stats.shapiro(female_AKS)
            normality_male = stats.shapiro(male_AKS)
            
            # Homogeneity of variances test
            variance_test = stats.levene(female_AKS, male_AKS)
            
            # Perform t-test assuming normality and homogeneity of variances
            if normality_female.pvalue > 0.05 and normality_male.pvalue > 0.05 and variance_test.pvalue > 0.05:
                # Use independent samples t-test
                t_statistic, p_value = stats.ttest_ind(female_AKS, male_AKS)
                plt.text(0.5, -0.2, f'p (t-test) = {p_value:.4f}\n', ha='center', va='center', transform=plt.gca().transAxes, fontsize=10, color='black')
                
            else:
                # Use Welch's t-test or Mann-Whitney U test as appropriate
                if variance_test.pvalue <= 0.05:
                    t_statistic, p_value = stats.ttest_ind(female_AKS, male_AKS, equal_var=False)
                    plt.text(0.5, -0.2, f'p (Welch\'s t-test) = {p_value:.4f}\n', ha='center', va='center', transform=plt.gca().transAxes, fontsize=10, color='black')
                    
                else:
                    U_statistic, p_value = stats.mannwhitneyu(female_AKS, male_AKS)
                    plt.text(0.5, -0.2, f'p (Mann-Whitney U test) = {p_value:.4f}\n', ha='center', va='center', transform=plt.gca().transAxes, fontsize=10, color='black')
                    
            
            significance_label = '*' if p_value < 0.05 else 'n.s.'
            
            
            # Define bracket position and text position
            x1, x2 = 0, 1  # Positions for 'MALE' and 'FEMALE' on x-axis
            y, h, col = summed_values_sex_recorded['AKS'].max()+0.1 - 0.05, 0.02, 'black'  # Adjust the y-position and height of the bracket
            
            # Draw the significance bracket
            plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, color=col)
            plt.text((x1 + x2) * 0.5, y + h, significance_label, ha='center', va='bottom', color=col, fontsize=10)
            
            # Display the plot
            plt.tight_layout()
            
            filename_box = f"{plot_path}boxplot_plot_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
            plt.savefig(filename_box)
            plt.close()

        #############################################################################
        ## if female only available. Only scatterplot will be used and no boxplot ##
        ############################################################################
        if len(female_data_length_recorded)>0 and len(male_data_length_recorded)==0:
            
            female_data_length_recorded.to_csv(f'{csv_path}FEMALE_info_{species_length_filter}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_.csv',index=False)
            # Perform linear regression for female data
            slope_female, intercept_female, r_value_female, p_value_female, std_err_female = linregress(female_data_length_recorded['AKS'], female_data_length_recorded['Length'])
            
            # Plot scatter plot for female and male data
            plt.figure(dpi=300)
            plt.scatter(female_data_length_recorded['AKS'], female_data_length_recorded['Length'], edgecolor='black', facecolor='white', marker='D' , label='Female')
            
            # Plot trend line for female data
            plt.plot(female_data_length_recorded['AKS'], slope_female*female_data_length_recorded['AKS'] + intercept_female, color='lightgrey',linestyle='solid', label='Female')
            
            # Annotate the plot with the equation of the trend line and p-value for female data
            plt.annotate(f'Female: y = {slope_female:.2f}x + {intercept_female:.2f} Pe = {p_value_female:.4f}', 
                         xy=(0.5, (plt.gca().xaxis.label.get_position()[1]-0.32)), xycoords='axes fraction',
                         fontsize=10, ha='center', va='top',
                         bbox=dict(boxstyle='round', fc='white', alpha=0.5))
            
            plt.title('Influence of size on leader-follower preference')
            
            plt.xlabel('AKS')
            plt.ylabel('Size(cm)')
            # Place the legend outside of the plot on the right side
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1),frameon=True, edgecolor='black',)
            plt.tight_layout()
            filename_scatter = f"{plot_path}scatter_plot_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
            plt.savefig(filename_scatter)
            
            plt.close()


        ########################################################################
        ## if only male available. Only scatterplot will be used and no boxplot
        #######################################################################
        if len(female_data_length_recorded)==0 and len(male_data_length_recorded)>0:
            
            # Perform linear regression for male data
            slope_male, intercept_male, r_value_male, p_value_male, std_err_male = linregress(male_data_length_recorded['AKS'], male_data_length_recorded['Length'])
            
            # Plot scatter plot for female and male data
            plt.figure(dpi=300)
            plt.scatter(male_data_length_recorded['AKS'], male_data_length_recorded['Length'], edgecolor='black', facecolor='darkgrey', marker='s', label='Male')
            
            # Plot trend line for male data
            plt.plot(male_data_length_recorded['AKS'], slope_male*male_data_length_recorded['AKS'] + intercept_male, color='darkgrey',linestyle='solid', label='Male')
            
            # Annotate the plot with the equation of the trend line and p-value for male data
            plt.annotate(f'Male: y = {slope_male:.2f}x + {intercept_male:.2f} P = {p_value_male:.4f}', 
                         xy=(0.5, (plt.gca().xaxis.label.get_position()[1]-0.42)), xycoords='axes fraction',
                         fontsize=10, ha='center', va='top',
                         bbox=dict(boxstyle='round', fc='white', alpha=0.5))
            
            plt.title('Influence of size on leader-follower preference')
            
            plt.xlabel('AKS')
            plt.ylabel('Size(cm)')
            # Place the legend outside of the plot on the right side
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1),frameon=True, edgecolor='black',)
            
            
            plt.tight_layout()  # Adjust the layout to make room for the legend
            filename_box = f"{plot_path}boxplot_plot_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
            plt.savefig(filename_box)
            plt.close()

    
    ########################################################################################################################################################
    ########################################################################################################################################################
    ########################################################################################################################################################
    
    #check if there is any output from the multiprocessing. If not only the allIDs csv is printed with an overview of the individuals
    if len(final_LF_info)>0:
        
        LF_grouped = final_LF_info.groupby('Name').agg({'L': 'sum', 'F': 'sum'}).reset_index()
    
        final_LF_info.to_csv(f'{csv_path}LF_info_{species_length_filter}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_.csv',index=False)
    
        LF_grouped.to_csv(f'{csv_path}LF_grouped_{species_length_filter}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_.csv', index=False)
    
    
    ################
    # Graphml output 
    ################

    # Create an empty graph to fuse the networks
    G_fused = nx.DiGraph()

    # Define a dictionary to store edge duplication frequency
    edge_frequency = {}  
    
    # Function to merge or update edge attributes
    def merge_edge_attributes(existing_attrs, new_attrs):
        # For each attribute in the new edge, add or update in the existing edge attributes
        for attr, value in new_attrs.items():
            if attr in existing_attrs:
                # Here you can define how to merge attributes if they already exist
                # For example, summing the values if they are numerical
                if isinstance(value, (int, float)) and isinstance(existing_attrs[attr], (int, float)):
                    existing_attrs[attr] += value
                else:
                    # If the attribute is not numerical, or you wish to handle it differently, adjust here
                    existing_attrs[attr] = value  # This line replaces the existing attribute
            else:
                # If the attribute does not exist, simply add it
                existing_attrs[attr] = value
        return existing_attrs
    
    # Iterate over all networks in the list and merge them into the fused graph
    for G in list_of_graphs:
        for u, v, data in G.edges(data=True):
            if G_fused.has_edge(u, v):
                # Merge edge attributes if the edge already exists
                existing_attrs = G_fused[u][v]
                merged_attrs = merge_edge_attributes(existing_attrs, data)
                G_fused[u][v].update(merged_attrs)
                # Increment frequency
                edge_frequency[(u, v)] = edge_frequency.get((u, v), 0) + 1
            else:
                # Add the edge with its attributes if it does not exist
                G_fused.add_edge(u, v, **data)
                edge_frequency[(u, v)] = 1
    
    # Set node attributes for the fused graph
    for node in G_fused.nodes():
        # Merge node attributes from individual graphs
        sex_attribute = {}
        size_attribute = {}
        for G in list_of_graphs:
            if node in G.nodes():
                sex_attribute[G.nodes[node]['sex']] = sex_attribute.get(G.nodes[node]['sex'], 0) + 1
                size_attribute[G.nodes[node]['measurement']] = size_attribute.get(G.nodes[node]['measurement'], 0) + 1
        # Set node attributes for the fused graph
        G_fused.nodes[node]['sex'] = max(sex_attribute, key=sex_attribute.get)
        G_fused.nodes[node]['measurement'] = max(size_attribute, key=size_attribute.get)
    
    # Set edge attributes for the fused graph
    for edge, freq in edge_frequency.items():
        G_fused.edges[edge]['weight'] = freq


    # Save the combined graph and save  as output
    nx.write_graphml(G_fused, f'{plot_path}combined_graph{species_length_filter}.graphml')

#########
# DAG TEST 
##########
        
    ## Function to get the unique edges in a directed graph. This ensures the correct counting of removed edges nececssary to create a DAG
    def get_unique_directed_edges(cycles):
        unique_edges = set()
        for cycle in cycles:
            # Create directed edges from the cycle
            edges = [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(len(cycle))]
            # Add to set to remove duplicates 
            unique_edges.update(edges)
        return list(unique_edges)
    
    # Function to remove edges to achieve DAG. It has the max_edges_to_remove parameter that will be used to check how many edges need to be removed to achieve DAG
    # The aim is always to have the lowest number possible.
    def remove_edges_to_achieve_dag(G, cycles, max_edges_to_remove):
        results = []
        removed_edges_count = 0

        # Check combinations of different sizes, up to the maximum specified
        for r in range(1, max_edges_to_remove + 1):
            for edges_to_remove in combinations(cycles, r):
                # Store the edges with their attributes before removing them
                edges_with_attrs = [(u, v, G[u][v]) for u, v in edges_to_remove]

                # Remove the edges
                G.remove_edges_from(edges_to_remove)
                
                # Check if the graph is now a DAG
                if nx.is_directed_acyclic_graph(G):
                    results.append((edges_to_remove, "makes G a DAG"))
                    removed_edges_count += 1
                # Add the edges back with their original attributes if not a DAG
                for u, v, attrs in edges_with_attrs:
                    G.add_edge(u, v, **attrs)

        if not results:
            results.append("No combination of up to {} edges removal can make G a DAG".format(max_edges_to_remove))
        
        #return the results sentence and the number of edges that were removed to acheive this. The latter will be used. 
        return results, removed_edges_count
    
    
    ############ part for the true sample ##############
    G_DAG_truesample =  copy.deepcopy(G_fused)
    G_null_model = copy.deepcopy(G_fused)
    
    #True sample count of removed edges to achieve DAG
    #start at 0 to have an ouput for a DAG in the raw graph
    DAG_counter = 0
    DAG_removed_edges_true_sample= ()

    #check if the graph is already a DAG
    if nx.is_directed_acyclic_graph(G_DAG_truesample):
        #if yes the removed edges will be set to the DAG counter which is 0
        DAG_removed_edges_true_sample = DAG_counter
        #DAG_counter will be set to number of edges to avoid the running of the following while loop 
        DAG_counter = G_DAG_truesample.number_of_edges()

    else:
        #if not DAG all edges of cycles will be saved and only the unique directed edges will be saved
        #this is done to avoid any double counting of edges that might overlap in multiple circles
        cycles_raw = list(nx.simple_cycles(G_DAG_truesample))
        cycles = get_unique_directed_edges(cycles_raw)

    while DAG_counter < G_DAG_truesample.number_of_edges():
        #DAG_counter is used to check for all possiiblities to remove the lowest number of edges to achieve DAG
        DAG_result = remove_edges_to_achieve_dag(G_DAG_truesample, cycles, max_edges_to_remove=DAG_counter)
        
        #If the DAG result is higher than 0 the number of removed edges will be saved and the DAG counter set to a value that breaks the while loop
        #It does not matter if for example the removal of two edges has multiple possibilities to create a DAG. Only the number 2 will be saved.
        if DAG_result[1] > 0:
            DAG_removed_edges_true_sample = DAG_counter
            DAG_counter = G_DAG_truesample.number_of_edges()
            
        DAG_counter +=1

    
    ##########################
    # Randomisation of network 
    ###########################
    
    # Create random network through edge swapping and crreate null model #####

    def swap_edges(G, num_swaps):
        G_swap = G.copy()
        edges = list(G_swap.edges())
        n_edges = len(edges)
        
        for _ in range(num_swaps):
            # Select two random edges
            idx1, idx2 = np.random.randint(0, n_edges, 2)
            edge1, edge2 = edges[idx1], edges[idx2]
            
            # Attempt swap if new edges don't already exist
            if not G_swap.has_edge(edge1[0], edge2[1]) and not G_swap.has_edge(edge2[0], edge1[1]):
                if edge1[0] != edge2[1] and edge2[0] != edge1[1]:  # Avoid self-loops
                    # Store attributes
                    attr1 = G_swap[edge1[0]][edge1[1]]
                    attr2 = G_swap[edge2[0]][edge2[1]]
                    
                    # Remove old edges
                    G_swap.remove_edge(edge1[0], edge1[1])
                    G_swap.remove_edge(edge2[0], edge2[1])
                    
                    # Add new edges with preserved attributes
                    G_swap.add_edge(edge1[0], edge2[1], **attr1)
                    G_swap.add_edge(edge2[0], edge1[1], **attr2)
                    
                    # Update edges list
                    edges[idx1] = (edge1[0], edge2[1])
                    edges[idx2] = (edge2[0], edge1[1])
        
        return G_swap

    ## DAG Null Model Analysis
    # Parameters
    n_permutations = Nruns
    num_swaps = G_null_model.number_of_edges() * 500  

    # Store results
    null_model_results = []

    # Run null model analysis
    for i in range(n_permutations):
        print(i)
        # Create randomized network
        G_null = swap_edges(G_null_model, num_swaps)
        
        # Calculate DAG edges for null model
        DAG_counter = 0
        if nx.is_directed_acyclic_graph(G_null):
            null_model_results.append(0)
            print()
        else:
            cycles_raw = list(nx.simple_cycles(G_null))
            cycles = get_unique_directed_edges(cycles_raw)
            
            while DAG_counter < G_null.number_of_edges():
                DAG_result = remove_edges_to_achieve_dag(G_null, cycles, max_edges_to_remove=DAG_counter)
                if DAG_result[1] > 0:
                    null_model_results.append(DAG_counter)
                    break
                DAG_counter += 1

    # Calculate the mean of the null distribution
    mean_null = sum(null_model_results) / len(null_model_results)

    # Calculate p-value
    p_value = sum(abs(result - mean_null) >= abs(DAG_removed_edges_true_sample - mean_null) for result in null_model_results) / n_permutations


    # plot the distribution
    plt.figure(dpi=300)
    unique_values, counts = np.unique(null_model_results, return_counts=True)

    # Convert counts to probability
    total_samples = np.sum(counts)
    probabilities = counts / total_samples

    plt.plot(unique_values, probabilities, '-o', label='Computed value\n (null model)', color='black')
    plt.axvline(DAG_removed_edges_true_sample, color='r', linestyle='--', 
                label='Computed value\n (combined network)')

    plt.xlabel('Minimum number of edges removed to convert network to DAG')
    plt.ylabel('Probability')
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    plt.legend()
    plt.title(f'Distribution of edge removals needed for DAG conversion\n in null model (n={Nruns})')

    #create a filename to save the output
    filename_DAG = f"{plot_path}DAG_test_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    # #save the output
    plt.savefig(filename_DAG)
    # plt.show()
    plt.close()
    
    
    #####################
    #  Hierachy analysis    
    #####################
    
    #Take a deepcopy of the fused network for the hierachy analysis 
    fused_network = copy.deepcopy(G_fused)
    #The number of simulations is based on the predined number of runs for all null models
    num_simulations = Nruns
       
    #Calculate the size of the largest connected component in the fused network
    size_lar_con_comp = len(max(nx.weakly_connected_components(fused_network), key=len))
    #Calculate the size of the largest out component+1 in the fused network
    largest_out_component_size = max(len(nx.descendants(fused_network, node)) + 1 for node in fused_network.nodes())
    
    
    ####################################################
    # Random network generation with random in- and outdegrees 
    ####################################################
    
    # Get the number of nodes and edges for the analysis
    Nodes_rand = fused_network.number_of_nodes()  # Number of nodes #originally at 21. But three individuals could have not been connected to the other individuals due to the time when they were tagged
    Edges_rand = fused_network.number_of_edges()  # Number of edge #soriginally at 19. But three individuals could have not been connected to the other individuals due to the time when they were tagged
    
    
    # Lists to store the sizes of the largest connected components and largest out-components
    largest_cc_sizes_rand = []
    largest_out_component_sizes_rand = []
    
    # Generate random directed networks and compute the required sizes
    for _ in range(num_simulations):
        G_directed_rand = gnm_random_graph(Nodes_rand, Edges_rand, directed=True)
        G_undirected_rand = G_directed_rand.to_undirected()
        
        # Largest connected component size
        largest_cc_rand = max(len(c) for c in connected_components(G_undirected_rand))
        largest_cc_sizes_rand.append(largest_cc_rand)
        
        # Largest out-component size
        largest_out_component_rand = max(len(descendants(G_directed_rand, node)) + 1 for node in G_directed_rand.nodes)
        largest_out_component_sizes_rand.append(largest_out_component_rand)
    
    #get the ratios
    ratios_rand = [out / cc for out, cc in zip(largest_out_component_sizes_rand, largest_cc_sizes_rand)]
    xy_rand = np.vstack([largest_cc_sizes_rand, ratios_rand])
    
    # Create KDE object
    kde_rand = gaussian_kde(xy_rand)
    
    # Evaluate KDE at the point (X=Size of largest connected component, Y=Size of largest out-component + 1 / Size of largest connected component)
    point_density_rand = kde_rand.evaluate([size_lar_con_comp, largest_out_component_size/size_lar_con_comp])
    
    
    # Plot the points on the plane
    plt.figure(figsize=(8, 6), dpi=300)
    plt.rcParams.update({'font.size': 14})  
    
    z = gaussian_kde(xy_rand)(xy_rand)
    
    # Sort the points by density, so the densest points are plotted last
    idx_orig = z.argsort()
    x, y, z = np.array(largest_cc_sizes_rand)[idx_orig], np.array(ratios_rand)[idx_orig], z[idx_orig]
    
    # Create the density contour plot
    plt.tricontourf(x, y, z, levels=14, cmap="Greys", alpha=0.7)
    plt.colorbar(label='Density')
    plt.scatter([size_lar_con_comp], [largest_out_component_size/size_lar_con_comp], color='Black', marker='*', s=200, label = "Combined network")
    plt.text(size_lar_con_comp + 0.5, largest_out_component_size/size_lar_con_comp, f'({point_density_rand[0]:.3f})', color='black', fontsize=12, verticalalignment='center')
    plt.xlabel('Size of largest connected component')
    plt.ylabel('Size of largest out-component + 1 /\n Size of largest connected component')
    plt.title('Random in- and outdegree sequences')
    plt.legend()
    plt.tight_layout()
    plt.xlim(min(x)-1, max(x)+1) 
    plt.ylim(min(y)-min(y)*0.1, max(y)+max(y)*0.1)  # x-axis starts from 7, extends to your max x + buffer (adjust as needed)
    
    #create a filename to save the output
    filename_hierarchy_rand = f"{plot_path}Hierachy_random_test_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    # #save the output
    plt.savefig(filename_hierarchy_rand)
    plt.close()
    
    
    
    ########################################################################################################
    # Random network generation preserved in- and outdegress sequences of the combined network 
    ########################################################################################################
    
    
    # Sort nodes numerically/alphabetically (depending on your node labels)
    node_list = sorted(fused_network.nodes(), key=lambda n: int(n) if str(n).isdigit() else str(n))
    
    # Get outdegree and indegree arrays in node order
    outdegrees = [fused_network.out_degree(n) for n in node_list]
    indegrees = [fused_network.in_degree(n) for n in node_list]
    
    # Lists to store the sizes of the largest connected components and largest out-components
    largest_cc_sizes_inout = []
    largest_out_component_sizes_inout = []
    
    # Generate random directed networks and compute the required sizes
    for _ in range(num_simulations):
        # Generate a directed graph with the given degree sequences
        G_directed_inout = directed_configuration_model(indegrees, outdegrees, create_using=nx.DiGraph)
        G_directed_inout = nx.DiGraph(G_directed_inout)  # Remove parallel edges and self-loops
        G_directed_inout.remove_edges_from(nx.selfloop_edges(G_directed_inout))
        G_undirected_inout = G_directed_inout.to_undirected()
        
        # Largest connected component size
        largest_cc_inout = max(len(c) for c in connected_components(G_undirected_inout))
        largest_cc_sizes_inout.append(largest_cc_inout)
        
        # Largest out-component size
        largest_out_component_inout = max(len(descendants(G_directed_inout, node)) + 1 for node in G_directed_inout.nodes)
        largest_out_component_sizes_inout.append(largest_out_component_inout)
    
    ratios_inout = [out / cc for out, cc in zip(largest_out_component_sizes_inout, largest_cc_sizes_inout)]
    
    # After generating largest_cc_sizes and ratios from your simulations
    xy_inout = np.vstack([largest_cc_sizes_inout, ratios_inout])
    
    # Create KDE object
    kde_inout = gaussian_kde(xy_inout)
    
    # Evaluate KDE at the point (X=Size of largest connected component, Y=Size of largest out-component + 1 / Size of largest connected component)
    point_density = kde_inout.evaluate([size_lar_con_comp, largest_out_component_size/size_lar_con_comp])
    
    # Plot the points on the plane
    plt.figure(figsize=(8, 6), dpi=300)
    plt.rcParams.update({'font.size': 14})  # Adjust the value as needed
    
    # Calculate the point density
    z_inout = gaussian_kde(xy_inout)(xy_inout)
    
    # Sort the points by density, so the densest points are plotted last
    idx_inout = z_inout.argsort()
    x_inout, y_inout, z_inout = np.array(largest_cc_sizes_inout)[idx_inout], np.array(ratios_inout)[idx_inout], z_inout[idx_inout]
    
    # Create the density contour plot
    plt.tricontourf(x_inout, y_inout, z_inout, levels=14, cmap="Greys", alpha=0.7)
    plt.colorbar(label='Density')
    plt.scatter([9], [8/9], color='Black', marker='*', s=200, label = "Combined network")
    plt.text(9 + 0.5, 8/9, f'({point_density[0]:.3f})', color='black', fontsize=12, verticalalignment='center')
    plt.xlabel('Size of largest connected component')
    plt.ylabel('Size of largest out-component + 1 /\n Size of largest connected component')
    plt.title('Preserving the in- and outdegree sequences\n of the combined network')
    plt.legend()
    plt.tight_layout()
    plt.xlim(min(x_inout)-1, max(x_inout)) 
    plt.ylim(min(y_inout)-min(y_inout)*0.1, max(y_inout)+max(y_inout)*0.1)  
    
    #create a filename to save the output
    filename_hierarchy_inout = f"{plot_path}Hierachy_inoutdegrees_test_{species_length_filter}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    # #save the output
    plt.savefig(filename_hierarchy_inout)
    plt.close()




######################################### END OF SCRIPT ###########################################################