#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#### Load packages and functions
library(shiny)
library(shinyMatrix)
mat = matrix(c(c(0, 0.33, 0.33, 0.33), #Transition freqs, given a base is picked
               c(0.33, 0, 0.33, 0.33),
               c(0.33, 0.33, 0, 0.33),
               c(0.33, 0.33, 0.33, 0)),
             ncol = 4)
nt <- c("A", "C", "G", "T")
colnames(mat) <- nt
rownames(mat) <- nt
mat_bases <- matrix(c(0.001055, 0.000526, 0.00045, 0.001607),
                    ncol = 4)
colnames(mat_bases) <- nt
####----------------------------------------------------------------------####


#### UI
ui <- fluidPage(
  sliderInput(inputId = "S_tot",
              label = "Number of synthetic sequences to generate per simulation cycle",
              min = 1,
              max = 10000,
              value = 1),
  sliderInput(inputId = "N_Cycles",
              label = "Number of simulation cycles to perform",
              min = 1,
              max = 1000,
              value = 300),
  textInput(inputId = "seq",
            label = "Template DNA sequence",
            placeholder = "only canonical DNA bases A, C, G, T"),
  numericInput(inputId = "seed",
               label = "Random number seed",
               value = 1,
               min = 1),
  matrixInput(inputId = "base_probs",
              label = "Base mutation probabilities. Determined by dNTP proportions in reaction.",
              value = mat_bases,
              rows = list(names = FALSE)),
  matrixInput(inputId = "base_transition_freqs",
              label = "Base transition frequencies from <row> to <column> (given a base mutates). Determined by polymerase and reaction conditions.",
              value = mat) ### Will have to adapt this to a vector later on
)
####----------------------------------------------------------------------####


#### Server
server <- function (input, output) {
  
  ### Input validation
  
  
  
  
  
  # *Input() functions
  # *Output() functions
}
####----------------------------------------------------------------------####


#### App
shinyApp(ui = ui, server = server)