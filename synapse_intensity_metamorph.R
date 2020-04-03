library(tidyverse)
library(ggpubr)
#######################################################################################
#######################################################################################
# ORGANIZING DATA
# use this script to get intensity values from Metamorph summary data
# use this before the script to get density per dendrite length
#######################################################################################
#######################################################################################

# choose metamorph summary file saved as .csv, skip first 3 rows
meta_data <-read.csv(file.choose(), skip = 3)   

# split Image.Time column into multiple columns to isolate a unique identifier (Minute) to act as proxy for image name
meta_data_sep<- meta_data %>% separate(Image.Time,  into = c("DateHour", "Minute","Sec", "MsecYear"), sep = ":")

# filters out the image header at the start of each new image 
# by filtering away rows containing any entries in column X.2
meta_data_tidy <- meta_data_sep %>% filter(meta_data_sep$X.2 == "")

# rename empty columns to add in "image" (the numerical order of image analysis) 
# and "overlay" (channel from which measurements were made 
# (with one region, 72 measurements are made (6 summary measurements * 12 overlays)))
meta_data_tidy <- rename(meta_data_tidy, Image = X, Overlay = X.1)

# filter data for only the columns of interest
meta_data_tidy <- meta_data_tidy %>% select (Image.Name, Minute, Summary.Statistic, Area,
                                             Average.intensity, Total.intensity, Image, Overlay)

# label the order of overlays/channels measured, repeated 6 times for the 6 summary measurements
meta_data_tidy$Overlay = rep(c('G','R','B','GR_G','GR_R','RB_R','RB_B',
                               'GB_G','GB_B','GRB_G','GRB_R','GRB_B'),each=6)

# return run length of hour column. Each cell/image should have 72 rows of data/measurements
rle <- rle(meta_data_tidy$Minute)
rle_entry_lengths <- data.frame(unlist(rle$lengths))

# check if all of them are 72 rows long 
# (sometimes fails if minute value is not unique)...maybe could think of a better check for this
warning(rle_entry_lengths==72)

# if all TRUE then its safe to auto label each "Image #" 
# below labels each image with its file name, pulled from the analysis/order folder 
# (ordered by date modified to match order of image analysis)

filename <- file.choose()
directory <- dirname(filename)
details=file.info(list.files(path = directory, full.names=TRUE))
details = details[with(details, order(as.POSIXct(mtime))), ]
files <- rownames(details)
meta_data_tidy$Image <- sapply(rep(files,each=72), basename)

# unblind the data, need to put CORRECT NAME of treatment and CORRECT NUMBER of cells in each condition.
# new column "Treatment" is created from these inputs (72 rows per image * number of images analysed per condition)

meta_data_tidy$Treatment <-  c(rep('CTL',(72*10)), rep('L7',(72*7)), rep('MK',(72*6)))

#######################################################################################
#######################################################################################
# GRAPHING
# filter tidyed dataset for what you're interested in graphing, then make graphs
#######################################################################################
#######################################################################################

# csv files gets read in as factors because it is a mix of numbers and text
# it cant be accurately graphed until corrected to be "as.numeric"

meta_data_tidy <- within(meta_data_tidy, {
  Average.intensity <- as.numeric(as.character(Average.intensity))
  Total.intensity <- as.numeric(as.character(Total.intensity))
  Area <- as.numeric(as.character(Area))
})

#create functions for making graphs to avoid having to write a bunch of code each time

#function to make boxplots:
#see ggplot2 to see all available customizations 

boxplotting <- function() {
  list(geom_boxplot(), 
       geom_point(position = position_jitterdodge(.05)),
       theme_classic(), 
       xlab(""),
       scale_fill_manual(values=c("grey51", "magenta", "cyan" ,"cyan3")),
       scale_x_discrete(labels=c("GRB_R" = "Shank3", "GRB_B" = "VGLUT1", "R" = "Shank3", "B" = "VGLUT1")))
}

#function to make barxplots:
#see ggbarplot and ggpar to see all available customizations 
barplotting <- function(df, x, y) {
  ggbarplot(df, x ,y,add = c("mean_se", "point"),
            #order = c('GRB_R', 'GRB_B'),
            color = "Treatment",  
            position = position_dodge()) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_color_manual(values=c("grey51", "dodgerblue1","dodgerblue2")) +
    scale_x_discrete(labels=c("GRB_R" = "Shank3", "GRB_B" = "VGLUT1", "R" = "Shank3", "B" = "VGLUT1"))
}


# filter the data to graph what you want. 
# here i'm filtering for "average" in the Summary.Statistc column ONLY for GRB "synaptic" overlays

graph_avg_total <- meta_data_tidy %>% filter(Summary.Statistic == "Average", 
                                             Overlay %in% c("GRB_R", "GRB_B"))

#force the order of appearance of variables on graph
graph_avg_total$Overlay <- factor(graph_avg_total$Overlay, levels = c("GRB_R", "GRB_B"))


#plot as a boxplot or a barplot, whichever you prefer:

#plot the Avg. Total Intensity of GRB/synaptic values as BAR plot using barplotting function
graph_avg_total %>% barplotting ("Overlay", "Total.intensity") + ylab("Avg. Total Intensity")

#plot as BOX plot using boxplotting function
graph_avg_total %>% ggplot(aes(x = Overlay, y = Total.intensity, fill=Treatment)) + boxplotting() +  ylab("Avg. Total Intensity")


# as a second example:
# filter the data for only SHANK3 (red channel) "average" in the Summary.Statistc (extrasynaptic&synaptic Shank3)
graph_avg_total_R <- meta_data_tidy %>% filter(Summary.Statistic == "Average", 
                                               Overlay %in% c("R"))
# graph this new data
graph_avg_total_R %>% barplotting ("Overlay", "Total.intensity") + ylab("Avg. Total Intensity")


#######################################################################################
#######################################################################################
#OPTIONAL DATA QUALITY CONTROL
#######################################################################################
#######################################################################################

#filter the data to exclude bad images with "!=" meaning "not" followed by the image name

graph_DROPPED <- meta_data_tidy %>% filter(Summary.Statistic == "Average",
                                           !Image %in% c('ctl_8.tif','L_15.tif','L_17.tif'),
                                           Overlay %in% c("GRB_R", "GRB_B"))


# graph the new synaptic data as a bar plot
graph_DROPPED %>% barplotting ("Overlay", "Total.intensity") + ylab("Avg. Total Intensity")


# get the new Shank3 data after excluding images

graph_DROPPED_RED <- meta_data_tidy %>% filter(Summary.Statistic == "Average", 
                                               !Image %in% c('ctl_8.tif','L_15.tif','L_17.tif'), 
                                               Overlay %in% c("R"))

#graph the new Shank3 data as a bar plot
graph_DROPPED_RED %>% barplotting ("Overlay", "Total.intensity") + ylab("Avg. Total Intensity")



#######################################################################################
#######################################################################################

# Optional: filter out cells with low numbers of synapses, set to whatever threshold you want for "count"

#######################################################################################
#######################################################################################

# get data for number of synapses
Count <- meta_data_tidy %>% filter(Summary.Statistic == "Count", Overlay %in% c("GRB_R")) %>%
  select(Area,Image,Treatment)  %>% 
  rename(Count = Area)


# get counts of puncta in RED channel

Count_R <- meta_data_tidy %>% filter(Summary.Statistic == "Count", Overlay %in% c("R")) %>%
  select(Area,Image,Treatment)  %>% 
  rename(Count = Area)

# get counts of puncta in BLUE channel

Count_B <- meta_data_tidy %>% filter(Summary.Statistic == "Count", Overlay %in% c("B")) %>%
  select(Area,Image,Treatment)  %>% 
  rename(Count = Area)


#graph counts per condition
#synaptic
barplotting(Count, "Treatment", "Count") + ylab("Synapse Number")
# all red channel
barplotting(Count_R, "Treatment", "Count") + ylab("Red Puncta Number")
# all blue channel
barplotting(Count_B, "Treatment", "Count") + ylab("Blue Puncta Number")

#######################################################################################

#filter the data to drop cells with less than "x" synapses (input as Area, in this example the cut off is 5)

Count_cutoff <- meta_data_tidy %>% filter(Summary.Statistic == "Count", Area > 5, 
                                          Overlay %in% c("GRB_R"))

# look at which cells are below the cut off
Cells_below_cutoff <- meta_data_tidy %>% filter(Summary.Statistic == "Count", Area < 5, 
                                                Overlay %in% c("GRB_R"))

Cells_below_cutoff
# exclude cells below the cutoff by image name 


graph_DROPPED <- meta_data_tidy %>% filter(Summary.Statistic == "Average",
                                           !Image %in% c('ctl_8.tif','L_17.tif','L_17max.tif','M_21.tif'),
                                           Overlay %in% c("GRB_R", "GRB_B"))

graph_DROPPED %>% barplotting ("Overlay", "Total.intensity") + ylab("Avg. Total Intensity")

############################################################################################
############################################################################################
#normalize in R
############################################################################################
############################################################################################
# once you have cleaned up your data and dropped bad cells, you may want to normalize the CTLs to 1 
# to be able to combine the data from replicates of experiments
# the following will convert the CTL average to 1 and the conditions relative to the control

# filter the data for just the CTL images and just the red channel
CTL_R_avg <-  graph_DROPPED %>% filter(Treatment == "CTL", Overlay %in% c("GRB_R"))

# create a new data frame to divide all CTL and values by the average of the CTL 

NormData_RED <- graph_DROPPED  %>% filter(Overlay %in% c("GRB_R"))
NormData_RED$NormAvgAVG = NormData_RED$Average.intensity/(mean(CTL_R_avg[["Average.intensity"]]))
NormData_RED$NormAvgTOT = NormData_RED$Total.intensity/(mean(CTL_R_avg[["Total.intensity"]]))


# check that the average of the CTL values is actually 1

ctl_norm_test <-  NormData_RED %>% filter(Treatment == "CTL", Overlay %in% c("GRB_R"))
mean(ctl_norm_test[["NormAvgTOT"]])


# repeat the above steps for the blue channel

CTL_B_avg <-  graph_DROPPED %>% filter(Treatment == "CTL", Overlay %in% c("GRB_B"))

NormData_BLUE <- graph_DROPPED  %>% filter(Overlay %in% c("GRB_B"))
NormData_BLUE$NormAvgAVG = NormData_BLUE$Average.intensity/(mean(CTL_B_avg[["Average.intensity"]]))
NormData_BLUE$NormAvgTOT = NormData_BLUE$Total.intensity/(mean(CTL_B_avg[["Total.intensity"]]))


ctl_norm_test <-  NormData_BLUE %>% filter(Treatment == "CTL", Overlay %in% c("GRB_B"))
mean(ctl_norm_test[["NormAvgTOT"]])


# combine the data frames of RED and BLUE normalized data into 1 data frame
all_norm <- bind_rows(NormData_RED, NormData_BLUE)

# look at this new data, see that CTL values avg to 1. should look exactly the same as graph_DROPPED data
barplotting(all_norm, "Overlay", "NormAvgTOT") + ylab("Avg. Total Intensity")

# export and save this data in the folder you are currently analyzing to have for later
# name the file whatever you want, here it is saved as "allnorm.csv"

write_csv(all_norm, path = file.path(dirname(directory),"allnorm.csv"))

