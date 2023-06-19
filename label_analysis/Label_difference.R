library(tidyverse)
library(ggplot2)
#Script to check and visualize difference between size of call versus distance between the labels for bionano calls filtered into the mastersheet

#Read in Mastersheet
callsheet <- read_csv("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Mastersheets/Label_Mastersheet.csv") 

#Add column of diff_size- The difference between the size of the call and the distance between the start and end label
callsheet <- callsheet %>% mutate(diff_size = (callsheet$Size/(callsheet$End-callsheet$Start)))
#Add column of diff_labels- difference in labels between the start and end labels in BP
callsheet <- callsheet %>% mutate(diff_labels = (callsheet$End-callsheet$Start))
#Selects coulmns of interest for better reading of tables in Rstudio if reading in different mastersheet
callsheet <- callsheet %>% select(Sample_ID,Id,Chromosome,Start,End,Size,Type,diff_size,diff_labels)

#Creates callsheet for looking at insertions as difference in diff_size values affects readability of graphs
callsheet_ins <- callsheet %>% filter(Type == "insertion")

#Boxplot to visualize the variablity of label spacing
ggplot(data=callsheet_del,mapping = aes(x=callsheet$Type,y=callsheet$diff_size)) +
  geom_boxplot()

#Filtering based on size to allow for closer visualization of different size boundries
callsheet_del <- callsheet %>% filter(Type == "deletion" & Size < 3000000 & Size >1500000)

#Makes column ID go to row. Mastersheet, just to check if one section of ID's is causing patterns in plots
callsheet_del <- callsheet_del %>% rowid_to_column()

#Point plot of deletion X-Size of call from Bionano, Y-diff_labels and colored by the percentage of size to diff_labels for better visualization or results
ggplot(data=callsheet_del,mapping = aes(x=callsheet_del$Size,y=callsheet_del$diff_labels)) +
  geom_point(aes(color = diff_size), size = 2) +
  scale_color_gradient2(low = "yellow", mid = "darkblue", high = "red", midpoint = mean(callsheet_del$diff_size))

#Point plot of deletion X-Size of call from Bionano, Y-diff_labels and colored by the percentage rowid to check for any pattern issues with calls
ggplot(data=callsheet_del,mapping = aes(x=callsheet_del$Size,y=callsheet_del$diff_labels)) +
  geom_point(aes(color = rowid), size = 2) +
  scale_color_gradient2(low = "yellow", mid = "darkblue", high = "red", midpoint = median(callsheet_del$rowid))


#Point plot of insertion x-Size of call y-diff_labels
ggplot(data=callsheet_ins,mapping = aes(x=callsheet_ins$Size,y=callsheet_ins$diff_labels)) + geom_point()

#Boxplot for insertion diff_labels
ggplot(data=callsheet_ins,mapping = aes(x=callsheet_ins$Type,y=callsheet_ins$diff_labels)) +
  geom_boxplot()

#Filtering of insertion by diff_labels to deal with possible outliers
callsheet_ins <- callsheet_ins %>% filter(diff_labels < 100000)

mean(callsheet_ins$diff_labels)
sum(callsheet_ins$diff_labels < 100000)
#Boxplot for insertion diff_labels
ggplot(data=callsheet_ins,mapping = aes(x=callsheet_ins$Type,y=callsheet_ins$diff_labels)) +
  geom_boxplot()
