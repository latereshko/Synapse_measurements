library(tidyverse)
library(ggplot2)


#######################################################################################
#######################################################################################
# Get density of synapses per dendrite length
# use this script after getting the intensity and count data from Metamorph summary data
#######################################################################################
#######################################################################################

#choose Metamorph Regions/Data file saved as .csv
meta_regions <-read_csv(file.choose())

#filter data to keep only the data from "Color Combine" channel (where the ROIs for dendrite lengths were drawn)
#BE CAREFUL THAT THE ORIGINAL POLYGON IS NOT INCLUDED IN THE DENDRITIC LENGTHS 
#(nothing in the Color Combine rows should match the area value of the Green, Red, Blue rows)

meta_regions_rename <- filter(meta_regions, `Image Name` == "Color Combine") %>% 
  select(1,3, Distance) %>% 
  rename(Channel = `Image Name`)

#split Image Date and Time column into multiple columns to isolate a unique identifier to act as proxy for image name
meta_regions_label <-separate(meta_regions_rename, `Image Date and Time`, into = c( "Time"), sep = 6)

#get lengths needed to repeat image name
rle <- rle(meta_regions_label$Time)
rle_entry_lengths <- as.numeric(unlist(rle$lengths))

#get order and names of images analysed ... 
#if running the script immediately after the summary data the filename will be the same
REGfilename <- filename
REGdirectory <- directory
REGdetails=file.info(list.files(path = REGdirectory, full.names=TRUE))
REGdetails = REGdetails[with(REGdetails, order(as.POSIXct(mtime))), ]
REGfiles <- rownames(REGdetails)
REGfiles <- REGfiles %>% tibble() %>% 
  rownames_to_column() %>% 
  rename(Order_Analysed = rowname, Image = "." )


#shorten pathname to basename
REGfiles$Image <- (basename(REGfiles$Image))

# repeat image name to label each row
meta_regions_label$Image <- rep(REGfiles$Image,(rle_entry_lengths))

# sum all dendrite lengths measured per image
distances <- meta_regions_label %>% select (Image, Distance) 
distances$Distance <- as.numeric(distances$Distance)

totallengths <- distances %>%
  group_by(Image) %>% 
  summarise(Distance = sum(Distance))


#join files by shared Image name
ordered_tot_lengths <- dplyr::full_join(REGfiles, totallengths, by = 'Image')


#join files by shared Image name - add Count from intensity analysis 
ordered_tot_lengths <- dplyr::full_join(Count, ordered_tot_lengths, by = 'Image')

#calculate density, put in new column
ordered_tot_lengths$Density <- ordered_tot_lengths$Count/ordered_tot_lengths$Distance

#graph as a bar SEM dots
ordered_tot_lengths %>% barplotting ("Treatment", "Density") 


#######################################################################################
#OPTIONAL
#filter the data to exclude bad images with "!=" meaning "not" followed by the image name
#######################################################################################
dens_DROPPED <- ordered_tot_lengths %>% filter(!Image %in% c('ctl_8.tif','L_15.tif','L_17.tif'))

#graph as a bar SEM dots
dens_DROPPED %>% barplotting ("Treatment", "Density") 

# save density data in the folder being analysed
write_csv(dens_DROPPED, path = file.path(dirname(directory),"dens_DROPPED.csv"))


