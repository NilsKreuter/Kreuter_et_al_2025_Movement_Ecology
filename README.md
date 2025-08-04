# Leader - follower analysis for acoustic telemetry data

## Overview

This repository contains the collection of scripts and workflows used to implement the leader-follower analysis presented in the manuscript "**Inferring leader-follower dynamics in three shark species using acoustic telemetry data**" by Nils Kreuter, Juan Fernández-Gracia, Víctor M. Eguíluz and Ana M. M. Sequeira published in Movement Ecology (2025).

## Scripts in this repository

-   **Leader_follower_method.py** : Main script to run the analysis. It calls all dependencies automatically. Only script the user needs to work with. Contains parameters for the analysis that the user needs to adjust.
-   **Leader_follower_KS** : Folder containing the python package and and functionalities used by the Leader_follower_method.py .
    -   **functions.py** : Located in the **Leader_follower_KS folder**. Contains all functions provided by the python package and that are called by **Leader_follower_method.py**.

## Installation

The python package can be installed via:

``` python3
pip install git+https://github.com/NilsKreuter/Kreuter_et_al_2025_Movement_Ecology.git
```

## Data and structure

The script is based on the data structure of the [IMOS ATF database](https://animaltracking.aodn.org.au/detection). Any data used for this method needs to follow this structure and column naming scheme to allow for the script to run properly.

**IMPORTANT:** Measurement data at the [IMOS ATF database](https://animaltracking.aodn.org.au/detection) are not always uniform (see example in table below). The **Leader_follower_method.py** has some code embedded at the beginning of the script to transform all measurements, from the three mentioned species in the publication, to a pure numeric representation of "Total length" in cm (exp: "TOTAL LENGTH = 3355 mm" —\> 335.5). This might need adjustment to fit your data.

| tag_id | detection_datetime | animal_sex | measurement | station_name | receiver_deployment_longitude | receiver_deployment_latitude |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| 123456 | 2020-01-16T13:52:12Z | FEMALE | TOTAL LENGTH = 3355 mm | Receiver_Island1 | 108.42 | -4.2 |
| 456343 | 2020-01-16T14:10:19Z | FEMALE | FORK LENGTH = 2153 mm | Receiver_Island1 | 108.42 | -4.2 |
| 234233 | 2020-01-16T14:11:00Z | MALE | TOTAL LENGTH = 244 cm | Receiver_Island1 | 108.42 | -4.2 |
| ... | ... | ... | ... | ... | ... | ... |

**Fish ID**: \
The script assigns unique character IDs to the provided numeric IDs. Both IDs are preserved in the output but the character IDs are used for every output to ease the identification of individuals.

**Receiver ID**:\
Receiver are assigned unique IDs (Receiver_1, Receiver_2, ...etc.) to allow for identification. The receiver IDs are based on unique combinations of receiver_deployment_longitude and receiver_deployment_latitude found in the data set. This avoids problematic naming schemes where receivers at different positions were assigned the same name. \

## How to use the main leader_follower script

**Leader_follower_method.py** is utilising multi-core processing approach to enhance analysis speed with larger datasets. The process on how to execute a multi-core processing python script is different depending on the environment (IDE or terminal). This might require specific file configuration in your IDE (for example: running the script in an external system terminal window when using [Spyder IDE](https://www.spyder-ide.org)).

### Setting the parameters {#setting-the-parameters}

In the "User input" section of the python script, several parameters need to be set for the analysis:

-   `cpus_used` : (numeric) Set the number of CPU cores used for analysis. The script uses multiprocessing to reduce analysis time. (Pay attention to specific requirements from your IDE and available cores.)

-   `csv_path_raw_data` : (path) Set the path to your .csv data file.

-   `timezone_raw_data` : (character) Set the timezone in which the data was saved.

-   `timezone_location` : (character) Set the timezone in which the data was collected.

-   `species_length_filter` : (character) Mention the species (this variable will be used for the plotting of graph titles and other output files.

-   `species_length_unit` : (character) Specifically for IMOS ATF data. Mention the length units mentioned in the data. Either "cm", "mm" or "cm and mm". The script will take care of the letters in the length measurements of the IMOS ATF data. (Conversion measurements are only available for tiger sharks, grey reef sharks and blacktip reef sharks). The code needs to be updated in case you use different data (See [**Data and structure**])

-   `pmax` : (numeric) Set the maximum p value for which a leadership pattern is determined (see manuscript for more details).

-   `Nruns` : (numeric) Set the number of reshuffled null model datasets (see manuscript for more details).

-   `min_int` : (numeric) Set the number of interactions that a pair must have to be considered for analysis (see manuscript for more details).

-   `output_path` : (path) Set the output path where the graphs and other output data should be saved (Exp: "C:/Users/Example/output/") .

### Example input

``` python3
#Set the path to your .csv data file
csv_path_raw_data = "C:/Users/Example/data.csv"
timezone_raw_data = "UTC"
timezone_location = "Australia/Brisbane"
species_length_filter = "Tiger"
species_length_unit = "mm"
pmax = 0.05
Nruns = 20
min_int = 25
output_path = "C:/Users/Example/output/"
```

# Output

**Leader_follower_method.py** will create a main folder and sub folder at the directory chosen for output path (see [**Setting the parameters**](#Setting-the-parameters)). The main folder is named after the selected species and the time of running the script (**Example**: Tiger_20220112_103549) while the sub folder is called "plots". The main folder will contain csv files as outputs, whereas "plots" will contain all plots created and a graphml object containing the aggregated network created for the selected species (See the manuscript for more details).

-   [**Main folder**]{.underline}
    -   **all_ID\_**...csv : Information about each individual in the dataset (ID_numeric, ID_character, measurement (cm), sex).
    -   **interaction_output\_**...csv : Number of edges and nodes sorted for each receiver in the data set, species tested and the number of minimum interactions set by the user.
    -   **LF_grouped\_**...csv : Each individual's (Name = ID_character) total counts as Leader (L) and Follower (F) across all receivers.
    -   **LF_info\_**...csv : Each individual's (Name = ID_character) counts as Leader (L) and Follower (F) for each receiver individually.
    -   **LF_strength\_**...csv : Combined information about Individuals (ID= ID_character), sex, length (total length in cm), their leader-follower event, its strength, which individual was either followed or lead (corresponding ID = cpd_ID) and location.
    -   **receiver_information**...csv : Longitude and Latitude coordinates for each receiver and its name.
-   [**"plots" subfolder**]{.underline}
    -   **boxplot_plot\_**...png : Box plot output for influence of sex on leadership preference (see manuscript for more details).
    -   **combined_graph**...graphml : Graphml file that contains all edge-, node-, individual information of the aggregated network.
    -   **DAG_test**...png : Output of the directed acyclic graph analysis (see manuscript for more details).
    -   **Hierachy_inoutdegrees**\_...png : Hierarchy analysis graph based on the preserved in- and outdegree sequences of the combined network (see manuscript for more details).
    -   **Hierachy_random**\_...png : Hierarchy analysis graph based on random in- and outdegree sequences of the combined network (see manuscript for more details).
    -   **detection_plot**...png : Overview of all detections for each individual in the entire dataset.
    -   **Network_plot_Receiver**...png: Leader-follower network for each receiver separately.
    -   **Receiver_1\_**...**AKS**\_plot.png: Output of the null model test for leader-follower dynamics for a pair. Black solid line represents the obtained AKS\* values of the null model and the red dashed line shows the AKS value from the real data (see manuscript for more details).
    -   **Receiver_1\_**...**TIME**\_plot.png: Overview of detections for each pair at the respective receiver. Includes number of interactions counted.
    -   **scatter_plot**...png : Scatter plot for sex-specific and overall influence of size on leadership preference (see manuscript for more details).

### Adjusting the node and edge size in the network plots

The size of nodes and edges might need to be adjusted depended on the number of individuals found to exhibit leader-follower behaviours at a single receiver.

Adjust the node size in the main leader-follower script with this parameter. \
(Lower number = smaller node)

``` python3
#create node sizes for the network
node_sizes = 1000
```

Adjust the edge size in the main leader-follower script by adjusting the numeric value (in this example **35**). \
(Lower number = smaller width of the edge)

``` python3
# # Draw edges and increasing the node size to reduce overlapping from nodes and edges
            nx.draw_networkx_edges(g, pos, edge_color="grey", arrows=True, arrowsize=10, node_size=node_sizes*2, 
                                   width=[data.get('A_KS', 1)*35 for _, _, data in g.edges(data=True)])
```

# Acknowledgement

Data was sourced from Australia’s Integrated Marine Observing System (IMOS) Animal Tracking Database ([https://animaltracking.aodn.org.au).](https://animaltracking.aodn.org.au).) IMOS is enabled by the National Collaborative Research Infrastructure Strategy (NCRIS). It is operated by a consortium of institutions as an unincorporated joint venture, with the University of Tasmania as Lead Agent. We acknowledge support from an ARC DP210103091 (awarded to A.M.M.S. and V.M.E.). The "Spangled Emperor at Heron Island" receiver array was made possible by fellowship DE120102459 from the Australian Research Council awarded to Alastair Harborne. We would like to thank Alastair Harborne and Michelle Heupel for their valuable feedback on the manuscript and thank all members of the “IMOS-ATF One Tree Island", "IMOS-ATF Heron Island" and "Spangled Emperor at Heron Island" projects for their work and their contribution to this study by making their data publicly available.

