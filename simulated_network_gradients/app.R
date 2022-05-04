
library(shiny)
library(tidyverse)
library(bslib)
library(patchwork)
source(here::here("_posts/2022-04-13-beta-diversity-and-species-responses/simulation_functions.R"))

ggplot2::theme_set(ggplot2::theme_bw() + theme(plot.title = element_text(size = 22)))

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty"),
    # Application title
    titlePanel("Simulating food webs"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          h3("Metaweb parameters"),
            numericInput("totalS", "Total species richness of Metaweb", min = 1, max = 70, 
                          step = 1, value = 25),
          h3("Species parameters"),
          
          h4("a, sensitivity (slope) of consumers to environment"),
          splitLayout(numericInput("avg_a",
                                   "Average",
                                   min = -6, max = 3, 
                                   step = .1, value = .3),
                      numericInput("sd_a", "SD", min = -3, max = 4, 
                                   step = .1, value = 1.1)
                      ),
          
          h4("b, Incidence of consumers"),
          splitLayout(numericInput("avg_b",
                                   "Average",
                                   min = -7, max = 3, 
                                   step = .1, value = -3),
                      numericInput("sd_b", "SD", min = 1, max = 4, 
                                   step = .1, value = .3)
                      ),
          
          h4("c, Degree distribution"),
          numericInput("sd_c", "SD of consumer degree distribution, logit scale", min = 1, max = 4, 
                       step = .1, value = .3),
          
          h4("correlations"),
          splitLayout(
            numericInput("corr_ab", "A & B", min = -1, max = 1, 
                         step = .1, value = .3),
            numericInput("corr_ac", "A & C", min = -1, max = 1, 
                         step = .1, value = .3),
            numericInput("corr_cb", "C & B", min = -1, max = 1, 
                         step = .1, value = .3)
            ),
          
          h3("sample information"),
          numericInput("nsites", "Number of bromeliads", min = 1, max = 100, 
                       step = 1, value = 20),
                            
             actionButton("simulate", "Simulate!")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           fluidRow(
             h2("Metaweb"),
             plotOutput("metaweb_figure"),
             ),
           fluidRow(
             h2("Bromeliad environments"),
             plotOutput("niche_figure"),
             plotOutput("richness_figure"),
             plotOutput("lcbd_figure"),
             plotOutput("Co_figure")
           )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  degrees <- eventReactive(input$simulate,{
    simulate_degrees(n_basal = 3, 
                     S = input$totalS,
                     degree_sd_logit = input$sd_p,
                     consumer_a_b_c = abc_values())
    })
  
  site_df <- eventReactive(input$simulate,{
    simulate_site_env(nsites = input$nsites)
    })
  
  abc_values <- eventReactive(input$simulate, {
    sim_abc(corr_ab = input$corr_ab,
            corr_ac = input$corr_ac,
            corr_cb = input$corr_cb, 
            sd_a = input$sd_a, sd_b = input$sd_b, sd_c = input$sd_c,
           S = input$totalS, n_basal = 3)
  })
  
  
  fake_data <- eventReactive(input$simulate, {
    req(site_df(), abc_values(), input$totalS, input$nsites)
    create_fake_data(nsites = input$nsites,
                     n_basal = 3, 
                     site_df = site_df(),
                     S = input$totalS,
                     consumer_a_b = abc_values(),
                     avg_a = input$avg_a,
                     avg_b = input$avg_b
                      )
  })
  
  simulated_network <- reactive({
    simulate_network(S = input$totalS, n_basal = 3, degrees = degrees())
  })
  

  simulated_observations <- reactive({
    simulate_presabs_obs(fake_data())
  })
  
  simulated_richness <- reactive({
    calculate_rich(simulated_observations())
  })
  
  simulated_sxs <- reactive({
    calculate_sxs(simulated_observations())
  })
  
  simulated_lcbd <- reactive({
    calculate_add_lcbd(sxs = simulated_sxs(), fake_rich = simulated_richness())
    
  })
  
  # Metaweb and bromeliad community combination
  matrix_subset_df <- reactive({
    make_matrix_subset_df(fake_data_obs = simulated_observations(),
                          n_basal = 3, 
                          one_network = simulated_network())
  })
  
  simulated_Co <- reactive({
    calculate_connect(matrix_subset_df = matrix_subset_df())
  })

# FIGURES -----------------------------------------------------------------

  
  
  output$metaweb_figure <- renderPlot({
    plot_degree_dist(degrees(), input$totalS) + 
      plot_network_mat(adj_mat = simulated_network(), S = input$totalS, n_basal = 3)
  })
  
  output$niche_figure <- renderPlot({
    
    site_df() |> 
      ggplot(aes(x = x, y = plogis(phi))) +  
      geom_point() + 
      labs(x = "Environment", y = "Average prob of species occuring there") + 
      
      fake_data() |> 
      ggplot(aes(x = env_x, y = plogis(logit_prob), group = sp)) + 
      geom_line() + 
      coord_cartesian(ylim = c(0,1)) + 
      labs(x = "Environment", y = "Probability of species occurrance") + 
      
      plot_layout(widths = c(1,2))
  })
  
  output$richness_figure <- renderPlot({
    simulated_richness() |> 
      ggplot(aes(x = env_x, y = richness)) + 
      geom_point() + 
      stat_smooth(method = "glm", method.args = list(family = "poisson")) + 
      labs(x = "Environment", y = "Species richness", title = "Richness is the sum of incidence")
  })

  output$lcbd_figure <- renderPlot({
    simulated_lcbd() |> 
      ggplot(aes(x = env_x, y = lcbd)) + 
      geom_point() + 
      stat_smooth(method = "gam", method.args = list(family = Gamma(link = "identity"))) +
      labs(x = "Environment", y = "LCBD", title = "LCBD can change with environment")
  })
  
  output$Co_figure <- renderPlot({
    simulated_Co() |> 
      ggplot(aes(x = env_x, y = Co)) + 
      geom_point() + 
      stat_smooth(method = "gam", method.args = list(family = "betar")) + 
      coord_cartesian(ylim = c(0,.4))+
      labs(y = "Connectance L/S^2", x = "Environment", title = "Connectance (L/S^2)")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
