library(shiny)
library(arrow)
library(DT)
library(shinythemes)
library(tidyverse)
library(grid)
library(GGally)
library(gridExtra)
library(feather)
library(jsonlite)
#library(plotly)

options(shiny.maxRequestSize = 30*1024^2)

load_tidy_feather <- function(path, generation='all'){
  df <- data.frame(arrow::read_feather(path))
  if (generation=='all'){
    df %>% 
      mutate(sumstat_name = paste0('sumstat_', sumstat_name),
             par_name = paste0('par_', par_name))  %>% 
      select(t, particle_id, distance, w,
             par_name, par_val, -index,
             par_name, par_val, -index,
             sumstat_name, sumstat_val) %>%
      spread(sumstat_name, sumstat_val) %>%
      spread(par_name, par_val) -> df
  }
  return(df)
}

return_names <- function(celltype, sim){
  # Return names for the summaryStatistic List given a celltype(str)
  # If Sim == True, Sim is appended else it's Exp
  if(sim == TRUE){
    mot <- paste0("sumstat_sumMot", celltype, "Sim")
    msd <- paste0("sumstat_sumMSD", celltype, "Sim")
  } else {
    mot <- paste0("sumstat_sumMot", celltype, "Exp")
    msd <- paste0("sumstat_sumMSD", celltype, "Exp")
  }
  return(c(msd, mot))
} 

combine_mot <- function(exp, sim, particle){
  sim <- sim %>% filter(particle_id == particle)
  
  target_Mot <- return_names("target", sim = TRUE)[[2]]
  
  sim_target   <- sim[[target_Mot]] %>% fromJSON()
  sim <- bind_rows('target' = sim_target,
                   .id = 'type')
  
  df <- bind_rows('Exp' = exp, 'Sim' = sim, .id = 'data')
  return(df)
}

combine_msd <- function(exp, sim, particle){
  sim <- sim %>% filter(particle_id == particle)
  
  target_MSD <- return_names("target", sim = TRUE)[[1]]
  #infected_MSD <- return_names("infected", col = col, sim = TRUE)[[1]]
  
  sim_target   <- sim[[target_MSD]] %>% fromJSON()%>%
    mutate(time_min = i*0.5) %>% select(-i)
  #sim_infected <- sim[[infected_MSD]] %>% fromJSON() %>%
   # mutate(time_min = i*0.5) %>% select(-i)
  sim <- bind_rows('target' = sim_target,
                   .id = 'type')
  
  df <- bind_rows('Exp' = exp, 'Sim' = sim, .id = 'data')
  return(df)
}

compare_mot <- function(exp, sim, particle){
  df_mot <- combine_mot(exp = exp,
                        sim = sim, 
                        particle = particle)
  plot <- df_mot %>% 
    gather(-data, -type, key = 'sumstat', value = 'value') %>% 
    ggplot(aes(x = type, y = value, col = data)) +
    facet_wrap(~sumstat, scale = 'free_y', 
               strip.position = "left", 
               labeller = as_labeller(c(speed = "Mean velocity [µm/min]", 
                                        angles = "Mean turning angle [rad]",
                                        arrest="Arrest coefficient", 
                                        straight="Straightness"))) + 
    geom_boxplot() +
    geom_point(position=position_jitterdodge(.1)) +
    ylab("") + xlab("") +
    theme_bw(base_size = 22) + 
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank(),
          legend.position = "none",
          panel.spacing = unit(2, "lines"))
  return(plot)
}

compare_msd <- function(exp, sim, particle){
  df_msd <- combine_msd(exp = exp,
                        sim = sim, 
                        particle = particle)
  plot <- df_msd %>% 
    ggplot(aes(x=time_min, y=mean/1000, col = data)) +
    geom_ribbon(aes(ymin=lower/1000,ymax=upper/1000, fill=data),alpha=0.3) +
    facet_wrap(~type, scales = "free_y", ncol = 1) + 
    geom_line() + theme_bw(base_size = 22) + 
    ylab("MSD [µm² x 10³]") + xlab("Time [min]") +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank(),
          panel.spacing = unit(3, "lines"),
          strip.text = element_text(size=22))
  return(plot)
}


my_density <- function(data, mapping, ...){
  # Using default ggplot density function
  # Using default ggplot density function
  ggplot(data = data, mapping = mapping) + 
    stat_density2d(aes(fill = ..density..), geom="tile", contour = FALSE) +
    scale_fill_gradientn(colours = viridis::viridis(100, option = "viridis")) +
    theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
}

my_point <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_point(shape = 21, alpha = 0.75, size = 0.5) 
}

matrix_plot <- function(df, t=NA){
  df <- df %>% select(-distance, -w, -particle_id, -starts_with('sumstat_'))
  if(!is.na(t)){
    df <- df %>% 
      filter(t==t) %>%
      select(-t)
  }
  ggpairs(df,
          lower = list(continuous = my_density),
          upper = list(continuous = my_point)) +
    theme_minimal(base_size = 22) +
    
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

compare_mot_msd <- function(exp_mot, exp_msd, sim, particle){
  p_mot <- compare_mot(exp_mot, sim, particle)
  p_msd <- compare_msd(exp_msd, sim, particle)
  grid.arrange(p_mot, p_msd, ncol = 2)
}

# Read data ####
df_mot <- target_mot_g
df_msd <- target_msd


ui <- fluidPage(
  titlePanel('HIV-FMC analysis'),
  sidebarLayout(
    sidebarPanel(
      selectInput('condition', 'Which condition?', 
                  choices=c('dense', 'loose')),
      selectInput('export', 'Which export option was chosen?', 
                  choices=c('last', 'all'), selected='last'),
      fileInput('featherdb', 'Select the feather file',
                accept=c('feather','.feather')),
      uiOutput('conditional_slider'),
      uiOutput('particle_id')),
    
    
    mainPanel(
      tabsetPanel(
        tabPanel('Table', DT::dataTableOutput('table')),
        tabPanel('Motility', plotOutput("motPlot",
                                        height = "600px")),
        tabPanel('ParameterMap', plotOutput("parPlot",
                                            height = "800px"))
      )
    )
  )
)

server <- function(input, output, session){
  
  # Read File ####
  input_file <- eventReactive(input$featherdb, {
    if (is.null(input$featherdb)) {
      return("Null")
    }
    # Read the text in the uploaded file
    load_tidy_feather(input$featherdb$datapath, generation = input$export)
  })
  
  # Conditional slider ####
  output$conditional_slider <- renderUI({
    conditionalPanel("input.export == 'all'",
                     sliderInput('generation', "Generation: ", 
                                 min = 1, 
                                 max=max(input_file()$t),
                                 value=1,
                                 step=1)
    )
  })
  
  # Particle Input
  # output$particle_id <- renderUI({
  #   sliderInput('particle_id','Particle:',
  #                value = input_file() %>%
  #                     filter(distance == min(distance)) %>%
  #                     pull(particle_id),
  #               )
  # })
  output$particle_id <- renderUI({
    shinyWidgets::sliderTextInput(inputId = "particle_id", 
                                  label = "Particle ID:", 
                                  choices = input_file() %>% select(particle_id, distance) %>%
                                    arrange(distance) %>% pull(particle_id))
  })
  
  
  # Table ####
  output$table <- DT::renderDataTable({input_file() %>% arrange(distance) %>%
      {if(input$export=='all') filter(., t == input$generation) else .} %>% 
      select(-starts_with("sumstat"))}, options=list(autoWidth = TRUE, scrollX = TRUE))
  
  # motPlot ####
  output$motPlot <- renderPlot({
    compare_mot_msd(exp_mot = df_mot,
                    exp_msd = df_msd,
                    sim = input_file(),
                    #col = input$condition,
                    particle = input_file() %>% 
                      filter(particle_id == input$particle_id) %>%
                      pull(particle_id))
  })
  
  # Parameter Plot ####
  output$parPlot <- renderPlot({
    matrix_plot(df = input_file(),
                {if(input$export=='all'){input$generation} else NA})
  })
  
}

shinyApp(ui=ui, server=server)
