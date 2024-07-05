
generate_plots <- function(plot_directory, directory = "~/"){
  
  # venn diagrams
  source(file.path(directory, "6GenerateFinalPlots/venn_diagrams.R"))
  create_venn_diagrams(plot_directory=plot_directory, directory=directory)
  
  source(file.path(directory, "6GenerateFinalPlots/volcano_plot.R"))
  create_volcano_plot(plot_directory=plot_directory, directory=directory)
  
  source(file.path(directory, "6GenerateFinalPlots/degree_plots.R"))
  create_degree_plots(plot_directory=plot_directory, directory=directory)
  
  
  source(file.path(directory, "6GenerateFinalPlots/create_graphs.R"))
  create_graphs(plot_directory=plot_directory, directory=directory)
}


generate_plots(plot_directory ='Plots', directory=file.path(getwd(), 'ROSMAP'))
