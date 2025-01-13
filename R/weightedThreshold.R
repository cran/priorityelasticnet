utils::globalVariables(c("sensitivity", "specificity", "Prediction", "Actual.Data", "Freq", "score", "X1", "X2", "xval", "yval", 
                         "result"))

#' A Shiny App for Model Evaluation and Weighted Threshold Optimization
#'
#' This function starts a Shiny application that enables users to interactively adjust the threshold for binary 
#' classification and view related metrics, the confusion matrix, ROC curve, and PR curve. 
#' The app also includes a feature for calculating the optimal threshold using a weighted version of Youden's J-statistic.
#' @param object A result from priorityelasticnet function with binomial model family.
#' @param ... Additional arguments
#' @details To calculate the optimal threshold, a weighted version of Youden's J-statistic (Youden, 1950) is used. 
#' The optimal cutoff is the threshold that maximizes the distance from the identity (diagonal) line. 
#' The function optimizes the metric (w * sensitivity + (1 - w) * specificity), where 'w' is the 
#' weight parameter adjusted using the second slider. After selecting the desired value on the optimal threshold slider, 
#' the user must press the "Set" button to update the threshold slider with the calculated optimal value. 
#' Metrics will then be automatically recalculated based on the user's selection. 
#' This function adapted from 'Monahov, A. (2021). Model Evaluation with Weighted Threshold Optimization 
#' (and the “mewto” R package). Available at SSRN 3805911.'
#' @importFrom shiny shinyApp fluidPage h2 strong fluidRow column wellPanel sliderInput tableOutput p tags tabsetPanel tabPanel plotOutput actionButton observeEvent updateSliderInput reactive renderTable renderPlot
#' @importFrom caret confusionMatrix  
#' @importFrom magrittr %>%
#' @importFrom pROC coords roc
#' @importFrom dplyr mutate rename ungroup summarise group_by across
#' @importFrom PRROC pr.curve
#' @importFrom plotrix draw.circle
#' @importFrom tidyr spread
#' @importFrom tibble rownames_to_column
#' @importFrom graphics lines
#' @importFrom broom tidy
#' @importFrom cvms plot_confusion_matrix font
#' @importFrom ggplot2 ggplot aes geom_line labs annotate scale_colour_gradient2 geom_segment geom_point theme_bw theme element_blank
#'
#' @return 
#'  No return value. This function is used for side effects only, specifically to launch a Shiny application for model evaluation with weighted threshold optimization.
#'  The Shiny app provides an interactive interface to visualize model performance metrics and optimize thresholds for classification models based on user-defined criteria.
#' 
#' @export
#' 
weightedThreshold <- function(object, ...) {
  
  if(!("priorityelasticnet" %in% class(object))){
    
    stop("The object must be of class 'priorityelasticnet'.")
    
  }
  
  if(object$family != "binomial"){
  
    stop("The model type (family) must be 'binomial'.")
    
  }
  
  
  actuals = object$actuals
  actuals[actuals == 1] = "yes"
  actuals[actuals == 0] = "no"
  probabilities = object$pred[,1]
  
  rocfun <- function(d_orig, d_pred) {
    roc_v <- suppressMessages(roc(response=d_orig,
                                  predictor=d_pred))
    return(roc_v)
  }
  
  confmat <- function(d_orig, d_pred, threshold) {
    d_pred_bin <- ifelse(d_pred>=threshold,"yes","no")
    confusionmatrix<-confusionMatrix(as.factor(d_pred_bin),as.factor(d_orig),
                                     positive="yes",
                                     dnn = c("Prediction", "Actual Data"))
    return(confusionmatrix)
  }
  
  comb <- function(rocfun_output, wght) {
    combinations <- as.data.frame((coords(roc = rocfun_output, x = rocfun_output$thresholds)))
    combinations <- combinations %>%
      mutate(youden = (wght * sensitivity + (1 - wght) * specificity))
    return(combinations)
  }
  
  prfun <- function(d_orig, d_pred) {
    pr_v <- pr.curve(scores.class0 = d_pred[d_orig=="yes"], # these are the predictions for default=no
                     scores.class1 = d_pred[d_orig=="no"], # these are the predictions for default=yes
                     curve=TRUE)
    
    return(pr_v)
  }
  
  shinyApp(
    
    ui <- fluidPage(
      h2(strong("Model Evaluation with Weighted Threshold Optimization"), style = "font-size:15px;"),
      fluidRow(
        column(4,
               wellPanel(
                 sliderInput("threshold", "Threshold:",
                             min = 0, max = 1,
                             value = 0.5)
               ),
               p(strong("Optimal Threshold"), style = "font-size:15px;"),
               wellPanel(
                 sliderInput("weight", "Weight on true positive rate maximization:",
                             min = 0, max = 1,
                             value = 0.5)
               ),
               tableOutput("values4"),
               fluidRow(column(9, tableOutput("values3")),
                        column(2, actionButton("set", "Set")))
        ),
        column(3,
               p(strong(" Confusion Matrix"), style = "font-size:15px;"),
               tableOutput("values2"),
               tableOutput("values")
        ),
        column(5,
               tags$head(
                 tags$style(type='text/css',
                            ".nav-tabs {font-size: 10px} ")),
               tabsetPanel(tabPanel("ROC Curve",
                                    plotOutput("plot1"),
                                    column(12, tableOutput("values5"), align="center")),
                           tabPanel("PR Curve", plotOutput("plot_PR"))
               )
        )
      )
    ),
    
    server <- function(input, output, session) {
      
      confmatL <- reactive({ confmat(actuals, probabilities, input$threshold) })
      
      rocfunL <- reactive({ rocfun(actuals, probabilities) })
      
      prfunL <- reactive({ prfun(actuals, probabilities) })
      
      combL <- reactive({ comb(rocfunL(), input$weight) })
      
      sliderValues <- reactive({
        rownames_to_column(as.data.frame(confmatL()$byClass), "Indicator") %>% rename(Value = 2)
      })
      
      sliderValues2 <- reactive({
        options(dplyr.summarise.inform = FALSE)
        as.data.frame(confmatL()$table) %>%
          group_by(Prediction,Actual.Data) %>%
          summarise(score=mean(Freq)) %>%
          spread(Actual.Data,score) %>%
          ungroup %>%
          mutate(across(is.numeric, as.integer)) %>%
          rename(Pred_Act = 1)
      })
      
      sliderValues3 <- reactive({
        results <- data.frame(c("Optimal threshold"),
                              c(combL()$threshold[which.max(combL()$youden)]))
        results <- results %>% rename(Method = 1, Threshold = 2)
      })
      
      observeEvent(input$set,{
        updateSliderInput(session,'threshold',value = combL()$threshold[which.max(combL()$youden)])
      })
      
      sliderValues4 <- reactive({
        weighttb <- as.data.frame(t(data.frame(c("Specificity (minimize FPR)",1-input$weight,"-","Sensitivity (maximize TPR)",input$weight))))
        weighttb <- weighttb %>%
          mutate(across(is.numeric, as.integer))
      })
      
      sliderValues5 <- reactive({
        roc_auc <- data.frame(c("AUC"),
                              c(rocfunL()$auc))
      })
      
      output$values <- renderTable({
        sliderValues()
      })
      
      output$values2 <- renderTable({
        sliderValues2()
      })
      
      output$values3 <- renderTable({
        format(sliderValues3(), nsmall = 4)
      }, colnames = FALSE)
      
      output$values4 <- renderTable({
        sliderValues4()
      }, digits = 0, colnames = FALSE)
      
      output$values5 <- renderTable({
        sliderValues5()
      }, colnames = FALSE)
      
      output$plot1 <- renderPlot({
        TPR <- (confmatL()$table[1,1])/(confmatL()$table[1,1]+confmatL()$table[2,1])
        FPR <- 1 - (confmatL()$table[1,2])/(confmatL()$table[1,2]+confmatL()$table[2,2])
        plot(rocfunL(), legacy.axes = TRUE, xlab="FPR - False Positives Rate", ylab="TPR - True Positives Rate")
        lines(x=c(TPR,TPR), y=c(0,FPR), col="cornflowerblue",lwd=2)
        lines(x=c(1,TPR), y=c(FPR,FPR), col="cornflowerblue",lwd=2)
        draw.circle(TPR,FPR,0.01,border="cornflowerblue",col="cornflowerblue")
      })
      
      output$plot_PR <- renderPlot({
        prec <- confmatL()$byClass[5]
        reca <- confmatL()$byClass[6]
        circ <- data.frame(xval = reca, yval=prec)
        scal <- max(1,as.integer(nrow(prfunL()$curve)/3000))
        ggplot(data.frame(prfunL()$curve[seq(1, nrow(prfunL()$curve), scal),]),aes(x=X1,y=X2)) +
          geom_line(size=1) + labs(x="Recall",y="Precision") +
          annotate("text", x = .0, y = .05, label = sprintf("AUC-PR Logit = %0.2f", prfunL()$auc.integral),hjust=0, size=4,
                   color="black", fontface=1) +
          scale_colour_gradient2(low="yellow", mid="green",high="blue") +
          geom_segment(aes(y = 0, x=reca), yend = prec, xend = reca, color = "cornflowerblue")+
          geom_segment(aes(y = prec, x=0), yend = prec, xend = reca, color = "cornflowerblue")+
          geom_point(data=circ, mapping=aes(x = xval, y = yval, size=15), color="cornflowerblue")+
          theme_bw(base_size = 16)+
          theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      })
      
      output$plot_confusionmatrix <- renderPlot({
        cfm <- tidy(confmatL()$table)
        plot_confusion_matrix(cfm,
                              target_col = "Actual Data",
                              prediction_col = "Prediction",
                              counts_col = "n",
                              font_counts = font(
                                size = 5,
                                angle = 0,
                                color = "black"
                              ),
                              font_normalized = font(
                                size = 5,
                                angle = 0,
                                color = "black"
                              ),
                              rm_zero_text=FALSE,
                              palette = "Blues"
        )
      })
    }
  )
}

