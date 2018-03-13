#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Old Faithful Geyser Data"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
       sliderInput("bins",
                   "Number of bins:",
                   min = 1,
                   max = 50,
                   value = 30),
       checkboxInput("box1", "show/hide whatever", value=TRUE),
       submitButton("Submit"),
       textInput("input1", "Enter tab 1", value = "Tab 1!"),
       textOutput("slopeOut"),
       textOutput("intOut")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("distPlot", brush = brushOpts(id = "brush1")),
       h3("Geiser Data for kids"),
       textOutput("text3"),
       tabsetPanel(type = "tabs",
                   tabPanel("Tab A", br(), textOutput("out1")),
                   tabPanel("Tab B", br(), textOutput("out2")))
    )
  )
))
