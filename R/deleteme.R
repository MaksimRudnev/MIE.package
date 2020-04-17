library(shiny)
ui <- fluidPage(mainPanel(
  plotOutput("plot1", click="plot_click"),
  plotOutput("plot2")
))
server <- function(input, output) {
  
  dt<- reactiveValues(plot = NULL)
  
  output$plot1 <- renderPlot({
    g<-ggplot2::ggplot(d, aes(x,y))+ggplot2::geom_point()
    dt$plot <- g
    g
  })
  
  output$plot2 <- renderPlot({
   
    df <- nearPoints(dt$plot$data,input$plot_click,"x","y", maxpoints=1, allRows=F, threshold=5)
    print(df)
    
    dt$plot
  })
}
x <- rnorm(100)
y <- rnorm(100)
z <- rnorm(100)
d <- data.frame(x=x,y=y, z=z)
shinyApp(ui = ui, server = server)