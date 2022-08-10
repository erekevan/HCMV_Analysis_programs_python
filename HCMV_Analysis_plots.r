#! /usr/local/bin/Rscript
# Author: Ahmed Al Qaffas
# Program: create hitsogram of FRF 

# Install package and load it. 
load_me <-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  } # Function was originally writen by (https://stackoverflow.com/users/4125693/matthew) 
}
options(warn = -1)
suppressMessages(load_me('ggplot2', 'dplyr', 'htmlwidgets' ,'plotly', 'argparse', 'filesstrings'))

#Functions:
## Plotting Main:
  # Task: 
      # Create a plot of number of reads on X position
  # Inputs: Data-frame; table with two cols where the first is number of reads and the second is the position (required) 
  #       : plot Title (optional)
  # Outputs: Static plot (jpeg)
  #        : Interactive plot (html)
plot_me <- function(df, plot_title='', break_point=5, ticks = 5) { 
  fig <- ggplot(df,
         aes(x=Position, y=number_Of_Reads,
             color=ifelse(number_Of_Reads>=break_point, 'red', 'black')
             )) + 
    geom_point(size=2, shape=21) +
    scale_color_identity() + 
    geom_col() +   
    geom_hline(yintercept = 5, linetype="dashed", color = "blue") +
    scale_y_continuous(name="Number of Reads", limits=c(0, max(df$number_Of_Reads + 2)),
                       breaks = seq(0, max(df$number_Of_Reads + 2), ticks)) +
    scale_x_continuous(name = 'Genomic Position')+
    ggtitle(plot_title) + 
    geom_text(aes(label=ifelse(number_Of_Reads> 5, paste(number_Of_Reads,Position, sep = ', '),''), hjust= 'inward', vjust= -1),
              check_overlap = TRUE, color = 'black') +
    theme_classic()
  return(fig)
}

# Plotting dynamic 
## Create an interactive plot with plotly express, this image can show overhover info and zoom in and out. 
plot_dy <- function(myfile, title='add_title'){
  dfig <- plot_ly(type = 'scatter', mode = 'markers') 
  dfig <- dfig %>% 
    add_trace(
      y = myfile$number_Of_Reads,
      x = myfile$Position, 
      marker = list(color=ifelse(myfile$number_Of_Reads>=5, 'red', 'black')),
      hoverinfo = 'y',
      hovertemplate = paste('Position: %{x:.2f}',
                            '\nNumber of Reads: <b>%{y}</b><extra></extra>'),
      showlegend = F
    ) %>% 
    layout(title = title)
  return(dfig)
}

# Function to check if argument input was provided
check_me <- function(x, y){
  if (is.null(x)){
    return(y)
  }
  return(x)
}

  
######################################################## Running pipeline ######################################################## 
# Suppressing Warning Text
# Enable system argument as input 
parser <- ArgumentParser(description='Generate Frequncy plots of HCMV mapped reads.')
parser$add_argument('-c', '--coords', type="character", help='file with reads coords and frequncy')
parser$add_argument('-o', '--outputname', type="character", help='Save file name as str (default: output)')
parser$add_argument('-ca', '--c_right_terminal', type="integer", help='c sequnce starting position')
parser$add_argument('-ab', '--ab_left_terminal', type="integer", help='a sequnce ending position')
parser$add_argument('-js', '--junction_start', type="integer",
                    help='Junction region coords (start of b seq)')
parser$add_argument('-je', '--junction_end', type="integer",
                    help='Junction region coords (end of c seq)')
args <- parser$parse_args()
if (is.null(args$coords)){
  cat('Frequncies file was not povided!\n')
  cat('Now existing...!!!...\n')
  quit()
}

# Loading the data file
#myfile <- read.delim(as.character('dist.txt'), sep = '', col.names = c('number_Of_Reads','Position'))
#save_file_code <- unlist(strsplit('dist.txt' , split = '[.]'))[1]
myfile<-read.delim(args$coords, sep = '', col.names = c('number_Of_Reads','Position'))
save_file_code <- check_me(args$outputname, 'output')


# # Initialize coords of the regions of inquiry 
cat(' The program uses the following coords unless specfic coords was added to the argument
    : ab ending position (left) = 2000
    : ca starting postion (right) = 233500
    : junctions coords:
        starting from 194252
        ending with 197770\n')
ab_coord <- check_me(args$ab_left_terminal, 2000)
jun_start_coord <- check_me(args$junction_start, 194252)
jun_ending_coord <- check_me(args$junction_end, 197770)
ca_coords <- check_me(args$c_right_terminal, 233500 )
dir.create(paste0(save_file_code,'_plots'))
setwd(paste0(save_file_code,'_plots/'))
# Create Whole Genomic Sequence mapping plots:
wgs_static_plot <- plot_me(myfile,plot_title ="Whole-Genome Alignment") # Static plot 
ggsave(plot= wgs_static_plot, filename = paste0(save_file_code, '_wgs.jpg'), device = "jpg") 
fig_wgs <- plot_dy(myfile, 'Whole-Genome Alignment') # Dynamic plot
htmlwidgets::saveWidget(as_widget(fig_wgs), paste0(save_file_code,".html"), selfcontained = FALSE)


### Create sub-region mapping plots:
## Only Static plots are outputted since the WGS dynamic plot has 'zoom' function. 
# Left Region:
ab_left_end_subset <- subset(myfile, myfile$Position <= ab_coord)
ab_left_end_plot <- plot_me(ab_left_end_subset, plot_title = 'Left Terminal Region', ticks = 1) # Static plot
ggsave(plot= ab_left_end_plot, filename = paste0(save_file_code, '_LR.jpg'), device = "jpg") 

# Junction region: 
junction_subset <- subset(myfile, myfile$Position  >= jun_start_coord & myfile$Position <=jun_ending_coord)
junction_plot <- plot_me(junction_subset, plot_title = 'Junction Region', ticks = 1)
ggsave(plot= junction_plot, filename = paste0(save_file_code, '_jun.jpg'), device = "jpg") 

# Right Region:
ca_subset <- subset(myfile, myfile$Position > ca_coords- 100)
ca_subset_plot <- plot_me(ca_subset, plot_title = 'Right Terminal Region', ticks = 5) # Static plot
ggsave(plot= ca_subset_plot, filename = paste0(save_file_code, '_RR.jpg'), device = "jpg") 






