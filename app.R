#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
#

list.of.packages <- c('shiny', 'shinydashboard', 'foreach', 'dplyr', 'reshape2', 'ggplot2', 'equatiomatic', 'zeallot', 'plotly', 'scales', 'xlsx2dfs', 'shinybusy', 'tibble', 'readr', 'hydroGOF', 'stringr', 'BiocManager')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

list.of.packages <- c("pmp", "sva")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

list.of.packages <- c("WaveICA2.0")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) devtools::install_github(dengkuistat/WaveICA_2.0,host="https://api.github.com")


library(shiny)
library(shinydashboard)
library(foreach)
library(dplyr)
library(reshape2)
library(ggplot2)
library(equatiomatic)
library(zeallot)
library(plotly)
library(scales)
library(xlsx2dfs)
library(shinybusy)
library(tibble)
library(readr)
library(pmp)
library(hydroGOF)
library(stringr)
library(WaveICA2.0)
library(sva)

options(shiny.maxRequestSize=50000*1024^2)
memory.limit(size = NA)
memory.size(max = NA)

aligningbatch_RTcorrected <- function(IS_RT_deviation, selectedbatch, aligningbatch, breaks, trend_type)
{
    equation <- list()
    colnum <- which(parse_number(gsub("[^0-9A-Za-z///' ]","." , selectedbatch ,ignore.case = TRUE)) == parse_number(colnames(IS_RT_deviation)))
    RTCorrectedbatch <- NA
    ### Create piecewise function for each batch
    for (i in 1:length(breaks))
    {
        if (i == 1)
        {
            IS_subset <- subset(IS_RT_deviation, IS_RT_deviation[,1] >0 & IS_RT_deviation[,1] <= breaks[i])
            data_subset <- subset(aligningbatch, aligningbatch$rt>0 & aligningbatch$rt <= breaks[i])
            
            if(trend_type[i] == "none")
            {
                data_subset$adjustedRT <- data_subset$rt
                RTCorrectedbatch <- rbind(RTCorrectedbatch, merge(data_subset, aligningbatch, by=c("mz","rt")))
            }
            else if(trend_type[i] == "linear")
            {
                equation[[i]] <- lm(as.formula(paste(colnames(IS_subset)[colnum],"~",colnames(IS_subset)[1])), data = IS_subset)
                RTvalues <- data.frame(data_subset$rt)
                colnames(RTvalues)[1] <- colnames(IS_subset)[1]
                pred_values <- predict(equation[[i]], newdata=RTvalues)
                data_subset$adjustedRT<-data_subset$rt+pred_values
                RTCorrectedbatch <- rbind(RTCorrectedbatch, merge(data_subset, aligningbatch, by=c("mz","rt")))
            }
            else if(trend_type[i] == "nonlinear")
            {
                equation[[i]] <- glm(as.formula(paste(colnames(IS_subset)[colnum],"~",colnames(IS_subset)[1])), family = poisson(), data = abs(IS_subset))
                RTvalues <- data.frame(data_subset$rt)
                colnames(RTvalues)[1] <- colnames(IS_subset)[1]
                pred_values <- predict(equation[[i]], newdata=RTvalues, type = "response")
                if(IS_subset[1,colnum]>0)
                {
                    data_subset$adjustedRT<-data_subset$rt+pred_values
                }
                else
                {
                    data_subset$adjustedRT<-data_subset$rt-pred_values
                }
                RTCorrectedbatch <- rbind(RTCorrectedbatch, merge(data_subset, aligningbatch, by=c("mz","rt")))
            }
            
            next;
        }
        
        else
        {
            
            IS_subset <- IS_RT_deviation[(findInterval(c(breaks[i-1],breaks[i]),IS_RT_deviation[,1] ))[1]:(findInterval(c(breaks[i-1],breaks[i]),IS_RT_deviation[,1] ))[2],]
            data_subset <- subset(aligningbatch, aligningbatch$rt>breaks[i-1] & aligningbatch$rt <= breaks[i])
            
            if(trend_type[i] == "none")
            {
                data_subset$adjustedRT <- data_subset$rt
                RTCorrectedbatch <- rbind(RTCorrectedbatch, merge(data_subset, aligningbatch, by=c("mz","rt")))
            }
            else if(trend_type[i] == "linear")
            {
                equation[[i]] <- lm(as.formula(paste(colnames(IS_subset)[colnum],"~",colnames(IS_subset)[1])), data = IS_subset)
                RTvalues <- data.frame(data_subset$rt)
                colnames(RTvalues)[1] <- colnames(IS_subset)[1]
                pred_values <- predict(equation[[i]], newdata=RTvalues)
                data_subset$adjustedRT<-data_subset$rt+pred_values
                RTCorrectedbatch <- rbind(RTCorrectedbatch, merge(data_subset, aligningbatch, by=c("mz","rt")))
            }
            else if(trend_type[i] == "nonlinear")
            {
                equation[[i]] <- glm(as.formula(paste(colnames(IS_subset)[colnum],"~",colnames(IS_subset)[1])), family=poisson(), data = abs(IS_subset))
                RTvalues <- data.frame(data_subset$rt)
                colnames(RTvalues)[1] <- colnames(IS_subset)[1]
                pred_values <- predict(equation[[i]], newdata=RTvalues, type = "response")
                if(IS_subset[1,colnum]>0)
                {
                    data_subset$adjustedRT<-data_subset$rt+pred_values
                }
                else
                {
                    data_subset$adjustedRT<-data_subset$rt-pred_values
                } 
                RTCorrectedbatch <- rbind(RTCorrectedbatch, merge(data_subset, aligningbatch, by=c("mz","rt")))
            }
            
            next;
        }
        
        return(RTCorrectedbatch, equation)
    } 
    
    RTCorrectedbatch <- filter(RTCorrectedbatch, rowSums(is.na(RTCorrectedbatch)) != ncol(RTCorrectedbatch))
    return(list("RTCorrected_batch"=RTCorrectedbatch, "correction_equation"= equation))
}

mz_pairs <- function(batch1_mz_rt, batch2_mz_rt, ppmError, RTDeviation)
{
    mz_pairs_batch12 <- matrix(ncol = 4)
    
    for (i in 1:nrow(batch1_mz_rt))
    {
        ppm <- abs((batch1_mz_rt[i,1]-batch2_mz_rt[,1])/batch1_mz_rt[i,1]*1000000)
        
        ppm <- ifelse(ppm < ppmError, 1, 0)
        
        rt_change <- abs(batch1_mz_rt[i,2]-batch2_mz_rt[,2])
        
        rt_change <- ifelse(rt_change < RTDeviation, 1, 0)
        
        ppm_rt <- ppm + rt_change
        
        index <- which((ppm_rt) == 2)
        
        if(length(index)>0)
        {
            if (length(index)==1)
            {
                for(j in 1:length(index))
                {
                    mz_pairs_batch12 <- rbind(mz_pairs_batch12, cbind(batch1_mz_rt[i,1],batch1_mz_rt[i,2],batch2_mz_rt[index[j],1],batch2_mz_rt[index[j],2]))
                }
            }
            
            else if (length(index)>=1)
            {
                #(print(i))
                batch2_matches <- matrix(ncol = 2)
                for(j in 1:length(index))
                {
                    if(nrow(mz_pairs_batch12) == 0 || (batch1_mz_rt[index[j],1] %in% mz_pairs_batch12[,1] & batch1_mz_rt[index[j],2] %in% mz_pairs_batch12[,2]) == FALSE || (batch2_mz_rt[index[j],1] %in% mz_pairs_batch12[,3] & batch2_mz_rt[index[j],2] %in% mz_pairs_batch12[,4]) == FALSE)
                    {
                        batch2_matches <- rbind(batch2_matches,cbind(batch2_mz_rt[index[j],1],batch2_mz_rt[index[j],2]))
                    }
                    else next
                    
                    batch2_matches <- na.omit(batch2_matches)
                    
                }
                if (nrow(batch2_matches) == 1)
                {
                    mz_pairs_batch12 <- rbind(mz_pairs_batch12, cbind(batch1_mz_rt[i,1],batch1_mz_rt[i,2],batch2_matches[1,1],batch2_matches[1,2]))
                }
                else if (nrow(batch2_matches) > 1)
                {
                    rt_diff <- (batch1_mz_rt[i,2]-batch2_matches[,2])
                    
                    if (any(abs(rt_diff)<0.01))
                    {
                        minrt <- which((abs(rt_diff)<0.01)==TRUE)
                        if (length(minrt) == 1)
                        {
                            mz_pairs_batch12 <- rbind(mz_pairs_batch12, cbind(batch1_mz_rt[i,1],batch1_mz_rt[i,2],batch2_matches[minrt,1],batch2_matches[minrt,2])) 
                        }
                        else
                        {
                            ppm_diff = (batch1_mz_rt[i,1]-batch2_matches[minrt,1])
                            k = which(min(abs(ppm_diff))==abs(ppm_diff))
                            mz_pairs_batch12 <- rbind(mz_pairs_batch12, cbind(batch1_mz_rt[i,1],batch1_mz_rt[i,2],batch2_matches[k,1],batch2_matches[k,2]))
                        }
                        
                    }
                    
                    else
                    {
                        k = which(min(abs(rt_diff))==abs(rt_diff))
                        mz_pairs_batch12 <- rbind(mz_pairs_batch12, cbind(batch1_mz_rt[i,1],batch1_mz_rt[i,2],batch2_matches[k,1],batch2_matches[k,2]))
                        
                        
                    }
                    
                }
                
            }
            
        }
    }
    
    
    duplicated_features <- unique(mz_pairs_batch12[which(duplicated(mz_pairs_batch12[,3:4])==TRUE),3:4])
    
    if (length(duplicated_features) == 2)
    {
        duplicated_features <- as.matrix(t(unique(mz_pairs_batch12[which(duplicated(mz_pairs_batch12[,1:2])==TRUE),1:2])))
        duplicated_rows <- mz_pairs_batch12[(mz_pairs_batch12[,3] %in% duplicated_features[,1] & mz_pairs_batch12[,4] %in% duplicated_features[,2]),]
        
        mz_pairs_batch12 <- mz_pairs_batch12[!(mz_pairs_batch12[,3] %in% duplicated_features[,1] & mz_pairs_batch12[,4] %in% duplicated_features[,2]),]
        
        for (i in 1:nrow(duplicated_features))
        {
            dups <- duplicated_rows[duplicated_rows[,3] %in% duplicated_features[i,1] & duplicated_rows[,4] %in% duplicated_features[i,2],]
            rt_diff = abs(dups[,4]-dups[,2])
            k = which(min(rt_diff[rt_diff<=0.02])==rt_diff)
            if (length(k) == 1)
            {
                mz_pairs_batch12 <- rbind(mz_pairs_batch12, dups[k,])
            }
            else
            {
                ppm_diff = (dups[,1]-dups[,3])
                k = which(min(abs(ppm_diff))==abs(ppm_diff))
                mz_pairs_batch12 <- rbind(mz_pairs_batch12, dups[k,])
            }
        }
    }else if (length(duplicated_features) > 2) 
    {
        duplicated_rows <- mz_pairs_batch12[(mz_pairs_batch12[,3] %in% duplicated_features[,1] & mz_pairs_batch12[,4] %in% duplicated_features[,2]),]
        
        mz_pairs_batch12 <- mz_pairs_batch12[!(mz_pairs_batch12[,3] %in% duplicated_features[,1] & mz_pairs_batch12[,4] %in% duplicated_features[,2]),]
        
        for (i in 1:nrow(duplicated_features))
        {
            dups <- duplicated_rows[duplicated_rows[,3] %in% duplicated_features[i,1] & duplicated_rows[,4] %in% duplicated_features[i,2],]
            rt_diff = abs(dups[,4]-dups[,2])
            k = which(min(rt_diff[rt_diff<=0.02])==rt_diff)
            if (length(k) == 1)
            {
                mz_pairs_batch12 <- rbind(mz_pairs_batch12, dups[k,])
            }
            else
            {
                ppm_diff = (dups[,1]-dups[,3])
                k = which(min(abs(ppm_diff))==abs(ppm_diff))
                mz_pairs_batch12 <- rbind(mz_pairs_batch12, dups[k,])
            }
        }
    }
    
    
    duplicated_features2 <- unique(mz_pairs_batch12[which(duplicated(mz_pairs_batch12[,1:2])==TRUE),1:2])
    
    if (length(duplicated_features2) == 2)
    {
        duplicated_features2 <- as.matrix(t(unique(mz_pairs_batch12[which(duplicated(mz_pairs_batch12[,1:2])==TRUE),1:2])))
        duplicated_rows <- mz_pairs_batch12[(mz_pairs_batch12[,1] %in% duplicated_features2[,1] & mz_pairs_batch12[,2] %in% duplicated_features2[,2]),]
        
        mz_pairs_batch12 <- mz_pairs_batch12[!(mz_pairs_batch12[,1] %in% duplicated_features2[,1] & mz_pairs_batch12[,2] %in% duplicated_features2[,2]),]
        
        for (i in 1:nrow(duplicated_features2))
        {
            dups <- duplicated_rows[duplicated_rows[,1] %in% duplicated_features2[i,1] & duplicated_rows[,2] %in% duplicated_features2[i,2],]
            rt_diff = abs(dups[,4]-dups[,2])
            k = which(min(rt_diff[rt_diff<=0.02])==rt_diff)
            if (length(k) == 1)
            {
                mz_pairs_batch12 <- rbind(mz_pairs_batch12, dups[k,])
            }
            else
            {
                ppm_diff = (dups[,1]-dups[,3])
                k = which(min(abs(ppm_diff))==abs(ppm_diff))
                mz_pairs_batch12 <- rbind(mz_pairs_batch12, dups[k,])
            }
        }
    }else if (length(duplicated_features2) > 2) 
    {
        duplicated_rows <- mz_pairs_batch12[(mz_pairs_batch12[,1] %in% duplicated_features2[,1] & mz_pairs_batch12[,2] %in% duplicated_features2[,2]),]
        
        mz_pairs_batch12 <- mz_pairs_batch12[!(mz_pairs_batch12[,1] %in% duplicated_features2[,1] & mz_pairs_batch12[,2] %in% duplicated_features2[,2]),]
        
        for (i in 1:nrow(duplicated_features2))
        {
            dups <- duplicated_rows[duplicated_rows[,1] %in% duplicated_features2[i,1] & duplicated_rows[,2] %in% duplicated_features2[i,2],]
            rt_diff = abs(dups[,4]-dups[,2])
            k = which(min(rt_diff[rt_diff<=0.02])==rt_diff)
            if (length(k) == 1)
            {
                mz_pairs_batch12 <- rbind(mz_pairs_batch12, dups[k,])
            }
            else
            {
                ppm_diff = (dups[,1]-dups[,3])
                k = which(min(abs(ppm_diff))==abs(ppm_diff))
                mz_pairs_batch12 <- rbind(mz_pairs_batch12, dups[k,])
            }
        }
    }
    
    return(as.data.frame(mz_pairs_batch12))
}

# Define UI for application that draws a histogram
ui <- dashboardPage(

    # Application title
    dashboardHeader(
        titleWidth = 300,
        title = "MS-MBAP"),
    
    dashboardSidebar(
        width = 300,
        sidebarMenu(
            menuItem("RT Correction", tabName = "RTCorrection", icon = icon("chart-area")),
            menuItem("Feature Alignment", tabName = "Alignment", icon = icon("object-ungroup")),
            menuItem("Missing Data Imputation", tabName = "Imputation", icon = icon("braille")),
            menuItem("Batch Effect Correction", tabName = "BECorrection", icon = icon("th")),
            menuItem("Additional Information", tabName = "Questions", icon = icon("question"))
        )),
        
    dashboardBody(
        tags$style(HTML(".sidebar-menu li a { font-size: 20px; }")),
        tabItems(
            ### RT Correction and Alignment Tab
            tabItem(tabName = "RTCorrection",
                    box(width = 4, height = 600, title = h3("File Inputs"),
                        fileInput('IS_File', 'Input Reference Compounds RT Table', accept=c("csv","comma-separated-values",'.csv')),
                        fileInput('gradient_file', 'Input Gradient Information', accept=c("csv","comma-separated-values",'.csv')),
                        actionButton("ID_RefBatch", "Identify Reference Batch", style="color: #fff; background-color: #449e48;border-color: #000000"),
                        br(),
                        br(),
                        fileInput('data', 'Input Excel File of Data', accept=c('.xlsx')),
                        actionButton("SelectRTChanges", "Select Batches for RT Correction", style="color: #fff; background-color: #449e48;border-color: #000000"),
                        add_busy_spinner(spin = "fading-circle")
                    ),
                    box(width = 8, height = 600,
                        h3("Reference Batch:",textOutput("medianBatch", inline = T)), br(),
                        plotlyOutput("RT_deviation_Plot")
                    ),
                    box(width = 12,
                        title = h3("Select Batches that Need RT Correction"),
                        checkboxGroupInput(inputId = "BatchCheck", label = "Choices",
                                           choices =  NULL
                        ), 
                            div(style="display: inline-block;",textInput("breaks", "Identify Breaks in the Data:", value = "", placeholder = "2, 12, 18")),
                            div(style="display: inline-block; width: 100px;",HTML("<br>")),
                            div(style="display: inline-block;",textInput("trendtype", "Identify the type of trend in RT change: (none/linear/nonlinear)", value = "", placeholder = "linear, linear, nonlinear")),
                        br(),
                        br(),
                        actionButton("PerformRTCorrection", "Perform RT Correction for Selected Batches", style="color: #fff; background-color: #449e48;border-color: #000000"),
                        downloadButton('downloaddata', label = "Download")
                    )
                    
            ),
            ### Feature Alignemnt Tab
            tabItem(tabName = "Alignment",
                    fluidRow(box(width = 6,
                                 title = h3("Align Features Across the Batches"),
                                 div(style="display: inline-block;",textInput("RTShift", "Retention Time Deviation (ex. 0.2)", value = "", placeholder = "0.2")),
                                 div(style="display: inline-block; width: 100px;",HTML("<br>")),
                                 div(style="display: inline-block;",textInput("ppmError", "ppm Error (ex. 10)", value = "", placeholder = "10")),
                                 div(style="display: inline-block; width: 100px;",HTML("<br>")),
                                 br(),
                                 br(),
                                 actionButton("AlignFeatures", "Align Features", style="color: #fff; background-color: #449e48;border-color: #000000"),
                                 downloadButton('download_all_alignedfeatures', label = "Download Aligned Features"),
                                 add_busy_spinner(spin = "fading-circle")
                    ),
                    box(width = 6,
                        title = h3("Allow Missing Data:"),
                        div(style="display: inline-block;",textInput("numMissing", "In how many batches can a feature be missing?", value = "", placeholder = "2")),
                        br(),
                        br(),
                        actionButton("Subset_AlignedFeatures", "Subset Aligned Features", style="color: #fff; background-color: #449e48;border-color: #000000"),
                        downloadButton('download_subset_alignedfeatures_data', label = "Download Data of Selected of Aligned Features"),
                        add_busy_spinner(spin = "fading-circle")
                    )),
                    box(width = 3, height = 600,
                        title = h3("Summary of Aligned Features:"),
                        tableOutput("alignedFeatures_summary")
                    ),
                    box(width = 9, height = 600,
                        title = h3("RT Deviation of Aligned Features:"),
                        plotlyOutput("Feaures_RT_deviation_Plot")
                    )
            ),
            
            
            ### Missing Data Imputation Tab
            tabItem(tabName = "Imputation",
                    box(width = 8, height = "900px",
                        title = h3("Create Simulated Data:"),
                        div(style="display: inline-block;","Simulate the missingness of data from the entire feature space within the feature space of no missing data"),
                        br(),
                        br(),
                        actionButton("simulate_missingness", "Simulate Missingness", style="color: #fff; background-color: #449e48;border-color: #000000"),
                        add_busy_spinner(spin = "fading-circle"),
                        plotlyOutput("modeledMissingnessplots", height = "700px")
                        ),
                    box(width = 4,
                        title = h3("Impute Missing Data:"),
                        plotlyOutput("NoMissingTrue",  height = "300px"),
                        selectInput("impute", "Select Imputation Method:", choices = c("small value", "mean", "median", "bayesian PCA", "kNN", "Random Forest")),
                        actionButton("impute_data", "Impute Data", style="color: #fff; background-color: #449e48;border-color: #000000"),
                        downloadButton('download_imputed_data', label = "Download Imputed Data"),
                        add_busy_spinner(spin = "fading-circle"),
                        plotlyOutput("ImputedData_Plot", height = "350px")
                    )
                    
                    ),
            
            ### Batch Effect Correction Tab
            tabItem(tabName = "BECorrection",
                    box(width = 6,  title = h3("File Inputs"),
                        fileInput('additionalinfo_file', 'Input Additional Sample Info', accept=c("csv","comma-separated-values",'.csv')),
                        br(),
                        actionButton("CorrectBatchEffect", "Perform Batch Effect Correction", style="color: #fff; background-color: #449e48;border-color: #000000"),
                        add_busy_spinner(spin = "fading-circle")
                    ),
                    box(width = 6, height = 300,
                        plotlyOutput("ImputedData_Plot2", height = "300px")
                    ),
                    box(width = 6, height = 600, title = h3("WaveICA 2.0 Corrected Data"),
                        plotlyOutput("WaveICACorrected", height = "450px"),
                        downloadButton('download_waveica_data', label = "Download WaveICA Corrected Data")
                    ),
                    box(width = 6, height = 600, title = h3("ComBat Corrected Data"),
                        plotlyOutput("ComBatCorrected", height = "450px"),
                        downloadButton('download_combat_data', label = "Download ComBat Corrected Data")
                    )
                    ),
            
            ### Additional Info Tab
            tabItem(tabName = "Questions",
                    fluidPage(
                        h3("Main Goal"),
                        h4("The object of this tool is data aggretation of multi-batch untargeted mass spectrometry data."),
                        br(),
                        h3("Upload 4 files to complete the analysis (see sample files for formatting)"),
                        h4("1) .csv with the rention time designation for all internal standards and known compounds"),
                        h4("2) .csv of your gradient information"),
                        h4("3) .xlsx of your data with each batch as a seperate tab"),
                        h4("4) .csv with batch, sample type, and injection order designation for each injection"),
                        br(),
                        h3("Retention Time Correction"),
                        h4("Using internal standards and known compounds that elute along the length of the run, identify a reference batch. Using the reference batch, model and correct retention times for deviating batches."),
                        br(),
                        h3("Feature Alignment"),
                        h4("Using a ppm error threshold and retention time deviation allowance, match features to the reference batch and aggregate."),
                        br(),
                        h3("Missing Data Imputation"),
                        h4("Identify the best algorithm to impute missing values by running a simulation on a subset of data which has no missing data and apply it to your entire dataset."),
                        br(),
                        h3("Batch Effect Correction"),
                        h4("Between WaveICA 2.0 and ComBat, identify the best batch effect correction algorithm for your data.")
                    )
                    )
        )
    )
)




# Define server logic
server <- function(input, output, session) {

    InternalStandards_RT <- reactive({
        read.csv(input$IS_File$datapath, header = TRUE, sep = ",", row.names = 1)
        #read.csv("C:/Users/meh9jb/Desktop/Gates/MultiBatchPipeline/IS_RTs.csv", header = TRUE, sep = ",", row.names = 1)
    })
    
    gradient_info <- reactive({
        read.csv(input$gradient_file$datapath, header = TRUE, sep = ",")
        #read.csv("C:/Users/meh9jb/Desktop/Gates/MultiBatchPipeline/Gradient.csv", header = TRUE, sep = ",")
    })
    
    complete_data <- reactiveValues(data = NULL)
    
    observeEvent(input$data$datapath,{
        show_modal_spinner(text = "Uploading Data") # show the modal window
        
        complete_data$data <- xlsx2dfs(input$data$datapath, rowNames = FALSE, colNames = TRUE)
        #complete_data$data <- xlsx2dfs("C:/Users/meh9jb/Desktop/Gates/MultiBatchPipeline/Complete_Data.xlsx", rowNames = FALSE, colNames = TRUE)
        
        remove_modal_spinner() # remove it when done
    })
    
    refbatch_variables <- reactiveValues( refBatch_value = NULL, IS_RT_deviation = NULL )
    
    observeEvent(input$ID_RefBatch,
                 {
                     show_modal_spinner(text = "Identifying Reference Batch")
                     medianBatch_value <- reactive({
                         ### Identify median batch
                         #Get Names of all Batches
                         
                         #test <<- InternalStandards_RT()
                         
                         ListBatch <- colnames(InternalStandards_RT())
                         
                         #For each compound (row), identify it's median RT across all batches
                         median_RTs<-apply(InternalStandards_RT(),1,median)
                         
                         #For each compound, identify in each batch if that batch's RT matches the median RT
                         median_yesno <- foreach(i=1:nrow(InternalStandards_RT()), .combine=rbind) %do% ifelse(InternalStandards_RT()[i,]==median_RTs[i], "yes", "no")
                         
                         #For each batch (column), caluculate how many compound are at the median value
                         ListBatch <- cbind(ListBatch, foreach(i=1:ncol(median_yesno), .combine=c) %do% sum(median_yesno[,i] == "yes")) 
                         
                         #Select the batch with the most median compounds as your reference batch
                         medianBatch <- ListBatch[,1][as.numeric(ListBatch[,2])==max(as.numeric(ListBatch[,2]))][1]
                         
                         return(medianBatch)
                     })
                     
                     refbatch_variables$refBatch_value <- medianBatch_value()
                     output$medianBatch <- renderText(medianBatch_value())
                     
                     IS_RT_deviation_matrix <- reactive({
                         ### Generate Retention Deviation Matrix 
                         #Create new matrix with Orignal RT - Refernce Batch RT to calcualte RT deviation
                         IS_RT_deviation <- foreach(i=1:ncol(InternalStandards_RT()), .combine = cbind) %do%  (InternalStandards_RT()[,i]-InternalStandards_RT()[,colnames(InternalStandards_RT())==medianBatch_value()])*-1
                         
                         #Rename IS_RT_deviation column names with batch names
                         colnames(IS_RT_deviation) <- colnames(InternalStandards_RT())
                         
                         #Rename IS_RT_deviation row names with compound names
                         rownames(IS_RT_deviation) <- rownames(InternalStandards_RT())
                         
                         #Replace reference batch column (which will be full of 0 values) with original RTs
                         IS_RT_deviation[,colnames(IS_RT_deviation)==medianBatch_value()] <- InternalStandards_RT()[,colnames(InternalStandards_RT())==medianBatch_value()]
                         
                         #Shift reference batch column to the beginning of the matrix (as column 1)
                         IS_RT_deviation <- as.data.frame(IS_RT_deviation[,c(which(colnames(IS_RT_deviation)==medianBatch_value()),which(colnames(IS_RT_deviation)!=medianBatch_value()))])
                         
                         return(IS_RT_deviation)
                     })
                     refbatch_variables$IS_RT_deviation <- IS_RT_deviation_matrix()
                     
                     ### Plot retention deviation matrix and gradient information
                     output$RT_deviation_Plot <- renderPlotly({
                         IS_RT_deviation_long <- reactive({
                             # Reformat data table to long version 
                             melt(IS_RT_deviation_matrix(), id.vars = c(1))})
                         
                         fig <- plot_ly(width = 950, height = 500)
                         fig <- fig %>% add_markers(data = IS_RT_deviation_long(), x=IS_RT_deviation_long()[,1], y= ~value, color = ~variable, colors =  hue_pal()(ncol(IS_RT_deviation_matrix())))
                         
                         y2 <- list(
                             overlaying = "y",
                             side = "right",
                             title = "Percentage of Solution B"
                             #,range = list(0,100)
                         )
                         
                         y3 <- list(
                             overlaying = "y",
                             side = "right",
                             anchor = "free",
                             position = 1,
                             range = list(0, 1),
                             title = "Flow Rate"
                         )
                         
                         fig <- fig %>% add_lines(data = gradient_info(), x = gradient_info()[,1], y = gradient_info()[,3], yaxis = "y2", name = "Gradient Line", line = list(color = 'rgb(142,142,142)', width = 2,dash = 'dash'))
                         
                         fig <- fig %>% add_lines(data = gradient_info(), x = gradient_info()[,1], y = gradient_info()[,2], yaxis = "y3", name = "Flow Rate", line = list(color = 'rgb(165,42,42)', width = 2, dash = 'dash'))
                         
                         hline <- function(y = 0, color = "red") {
                             list(
                                 type = "line",
                                 x0 = 0,
                                 x1 = 0.9,
                                 xref = "paper",
                                 y0 = y,
                                 y1 = y,
                                 line = list(color = color, dash="dot")
                             )
                         }
                         
                         fig <- fig %>% layout(
                             shapes = list(hline(0.2), hline(-0.2)),
                             title = paste("Known Compounds Retention Deviation Compared to",medianBatch_value()), 
                             yaxis2 = y2,
                             yaxis3 = y3,
                             xaxis = list(title="Retention Time (min)", domain = c(0, 0.9)),
                             yaxis = list(title="Retention Deviation (min)"),
                             legend = list(x = 1.1),
                             margin=list(t = 50)
                         )
                         
                         fig
                     })
                     remove_modal_spinner()
                 }, ignoreInit = TRUE)
    
    
    observeEvent(input$SelectRTChanges,
                 {
                     ##Remove Later
                     test <<- complete_data$data
                     
                     updateCheckboxGroupInput(session, "BatchCheck",
                                              label = "Choices",
                                              choices = names(complete_data$data)[which(names(complete_data$data) != refbatch_variables$refBatch_value)],
                                              inline   = TRUE
                     )
                 }, ignoreInit = TRUE)
    
    modifed_batches <- reactiveValues(listModifiedBatches = NULL )
    
    
    observeEvent(input$PerformRTCorrection,{
        
        show_modal_spinner(text = "Adjusting Retention Times for Selected Batches")
        
        if(length(input$BatchCheck) > 0)
        {
            modifed_batches$listModifiedBatches <<- input$BatchCheck
            
            for(i in 1:length(input$BatchCheck))
            {
                x <- match(input$BatchCheck[i], names(complete_data$data))
                breaks <<- as.numeric(unlist(strsplit(input$breaks, split = ",")))
                trend_type <<- unlist(strsplit(gsub(" ", "", input$trendtype, fixed = TRUE), split = ","))
                
                aligning_info <- complete_data$data[[x]][,1:2]
                
                modified_info <- aligningbatch_RTcorrected(refbatch_variables$IS_RT_deviation, input$BatchCheck[i], aligning_info, breaks, trend_type)
                
                complete_data$data[[x]] <- cbind(modified_info$RTCorrected_batch, complete_data$data[[x]][,-1:-2])
            }
            ##Remove Later
            test_modified <<- complete_data$data   
        }
        
        remove_modal_spinner()
    }, ignoreInit = TRUE)
    
    aligned_features_info <<- NULL
    
    allBatches_mzRTinfo <- reactive({
        allBatch_mzRTinfo <<- NULL
        
        for (i in 1:length(complete_data$data))
        {
            if (any(names(complete_data$data[i]) == modifed_batches$listModifiedBatches) == TRUE)
            {
                allBatch_mzRTinfo[i] <<- list(complete_data$data[[i]][,1:3])
            }
            
            else
            {
                allBatch_mzRTinfo[i] <<- list(complete_data$data[[i]][,1:2])
            }
        }
        
        names(allBatch_mzRTinfo) <<- names(complete_data$data)
        
        return(allBatch_mzRTinfo)
    })
    
    output$downloaddata <- downloadHandler(
        filename = function() {paste('featurelist_withAdjustedRTs', Sys.Date(), '.xlsx', sep = '')},
        content = function(file) {
            showModal(modalDialog("Dowloading Data with Adjusted RTs", footer = NULL))
            on.exit(removeModal())
            dfs2xlsx(allBatches_mzRTinfo(), file, rowNames = FALSE)}
    )
    
    alignedFeatures <- reactiveValues( all_aligned_features = NULL, aligned_features_subset = NULL, interested_features_dataframe = NULL, nomissingdata = NULL)
    
    observeEvent(input$AlignFeatures, {
        
        start_time <- Sys.time()
        show_modal_progress_line(text = "Matching Batch 1")
        
        allBatches_mzRTinfo <- allBatches_mzRTinfo()
        
        refbatch_index <- which(names(allBatches_mzRTinfo) == refbatch_variables$refBatch_value)

        aligned_features_info <<- allBatches_mzRTinfo[[refbatch_index]]
        
        colnames(aligned_features_info) <<- c(paste0(refbatch_variables$refBatch_value,"mz", sep = ""), paste0(refbatch_variables$refBatch_value,"rt", sp = ""))
        

             n <- length(complete_data$data)

            for (i in 1:n){

            update_modal_progress(value = i/n, text = paste("Matching Batch", i))

            if(ncol(allBatches_mzRTinfo[[i]]) == 3){

                mz_pair_info <- mz_pairs(allBatches_mzRTinfo[[refbatch_index]], allBatches_mzRTinfo[[i]][,c(1,3)], as.double(input$ppmError), as.double(input$RTShift))

                colnames(allBatches_mzRTinfo[[i]]) <- c(paste0(names(allBatches_mzRTinfo)[i],"mz", sep = ""), paste0(names(allBatches_mzRTinfo)[i],"rt", sep = ""), paste0(names(allBatches_mzRTinfo)[i],"editedrt", sep = ""))

                colnames(mz_pair_info) <- c(paste0(refbatch_variables$refBatch_value,"mz", sep = ""), paste0(refbatch_variables$refBatch_value,"rt", sp = ""), paste0(names(allBatches_mzRTinfo)[i],"mz", sep = ""), paste0(names(allBatches_mzRTinfo)[i],"editedrt", sep = ""))

                mz_pair_info <- merge(allBatches_mzRTinfo[[i]],mz_pair_info, by=c(paste0(names(allBatches_mzRTinfo)[i],"mz", sep = ""), paste0(names(allBatches_mzRTinfo)[i],"editedrt", sep = "")))
                

                aligned_features_info <- merge(aligned_features_info,mz_pair_info, by=c(paste0(refbatch_variables$refBatch_value,"mz", sep = ""), paste0(refbatch_variables$refBatch_value,"rt", sep = "")), all = TRUE)

            }

            else if (names(allBatches_mzRTinfo[i]) == refbatch_variables$refBatch_value) next


            else {
                mz_pair_info <- mz_pairs(allBatches_mzRTinfo[[refbatch_index]], allBatches_mzRTinfo[[i]],  as.double(input$ppmError), as.double(input$RTShift))

                colnames(allBatches_mzRTinfo[[i]]) <- c(paste0(names(allBatches_mzRTinfo)[i],"mz", sep = ""), paste0(names(allBatches_mzRTinfo)[i],"rt", sep = ""))

                colnames(mz_pair_info) <- c(paste0(refbatch_variables$refBatch_value,"mz", sep = ""), paste0(refbatch_variables$refBatch_value,"rt", sep = ""), paste0(names(allBatches_mzRTinfo)[i],"mz", sep = ""), paste0(names(allBatches_mzRTinfo)[i],"rt", sep = ""))

                aligned_features_info <- merge(aligned_features_info,mz_pair_info, by=c(paste0(refbatch_variables$refBatch_value,"mz", sep = ""), paste0(refbatch_variables$refBatch_value,"rt", sep = "")), all = TRUE)
            }


            }
             remove_modal_progress()
             show_modal_spinner(text = "Aggregating Batches")
            
             aligned_features_info <- aligned_features_info[rowSums(is.na(aligned_features_info)) != ncol(aligned_features_info), ]
             
            ##remove later 
            aligned_features_info <<- aligned_features_info
            
            alignedFeatures$all_aligned_features <- aligned_features_info[rowSums(is.na(aligned_features_info[ , grepl( "mz" , names( aligned_features_info))])) < n-1,] 
            
            all_aligned_features <<- alignedFeatures$all_aligned_features
            
            summary_info <- data.frame(cbind(names(summary(as.factor(rowSums(is.na(aligned_features_info[ , grepl( "mz" , names(aligned_features_info))]))))), unclass(summary(as.factor(rowSums(is.na(aligned_features_info[ , grepl( "mz" , names(aligned_features_info))])))))), check.names = FALSE)
            
            colnames(summary_info) <- c("# of Batches Missing In", "# of Features")
            
            output$alignedFeatures_summary <- renderTable(summary_info)
            
            end_time <- Sys.time()

            howlong <<- end_time - start_time
            
            remove_modal_spinner()

    }, ignoreInit = TRUE)
    
    observeEvent(input$Subset_AlignedFeatures, {
        
        show_modal_spinner(text = "Subetting Features")
        
        aligned_feature_info_subset <<- alignedFeatures$all_aligned_features[rowSums(is.na(alignedFeatures$all_aligned_features[ , grepl( "mz" , names(alignedFeatures$all_aligned_features))])) <= as.integer(input$numMissing),] 
        
        alignedFeatures$aligned_features_subset <- aligned_feature_info_subset
        
        interested_features_dataframe <- aligned_feature_info_subset
        
        for (i in 1:ncol(refbatch_variables$IS_RT_deviation))
        {
            colnums <- which(parse_number(gsub("[^0-9A-Za-z///' ]","." , colnames(refbatch_variables$IS_RT_deviation)[i] ,ignore.case = TRUE)) == parse_number(gsub("[^0-9A-Za-z///' ]","." , colnames(aligned_feature_info_subset) ,ignore.case = TRUE)))
            
            batch_index <- which(parse_number(gsub("[^0-9A-Za-z///' ]","." , colnames(refbatch_variables$IS_RT_deviation)[i] ,ignore.case = TRUE)) == parse_number(gsub("[^0-9A-Za-z///' ]","." , names(complete_data$data) ,ignore.case = TRUE)))
            
            interested_features_dataframe <- merge(interested_features_dataframe, complete_data$data[batch_index], by.x = c(colnums[1],tail(colnums,1)), by.y = 1:2, all.x = TRUE)
            
        }
        
        interested_features_dataframe <- interested_features_dataframe[ , !grepl( "adjustedRT" , names(interested_features_dataframe))]
        interested_features_dataframe <- interested_features_dataframe[ , !grepl( "editedrt" , names(interested_features_dataframe))]
        interested_features_dataframe <<- interested_features_dataframe
        alignedFeatures$interested_features_dataframe <- interested_features_dataframe
        
        output$Feaures_RT_deviation_Plot <- renderPlotly({
        
        feature_coord_info <- melt(aligned_feature_info_subset[ , grepl("rt", names(aligned_feature_info_subset)) & !grepl("edited", names(aligned_feature_info_subset))], id.vars = c(1))
        
        fig <- plot_ly(#width = 950, 
            height = 500)
        fig <- fig %>% add_markers(data = feature_coord_info, feature_coord_info[,1], y= ~(feature_coord_info[,1]-value), color = ~variable, colors =  hue_pal()(ncol(refbatch_variables$IS_RT_deviation)))
        
        fig <- fig %>% layout(
            title = paste("Retention Deviation of Aligned Features Compared to",refbatch_variables$refBatch_value), 
            xaxis = list(title="Retention Time (min)"),
            yaxis = list(title="Retention Deviation (min)"),
            legend = list(x = 1.1),
            margin=list(t = 50)
        )
        
        fig
    })
        remove_modal_spinner()
        }, ignoreInit = TRUE)
    
    output$download_all_alignedfeatures <- downloadHandler(
        filename = function() {paste("allalignedFeatures-", Sys.Date(), ".csv", sep = "")},
        content = function(file) {
            showModal(modalDialog("Dowloading All Aligned Features", footer = NULL))
            on.exit(removeModal())
            write.csv(alignedFeatures$all_aligned_features, file, row.names = FALSE)}
    )
    
    output$download_subset_alignedfeatures_data <- downloadHandler(
        filename = function() {paste("selected_alignedFeatures_data-", Sys.Date(), ".csv", sep = "")},
        content = function(file) {
            showModal(modalDialog("Dowloading Aggregated Data of Interested Features", footer = NULL))
            on.exit(removeModal())
            write.csv(alignedFeatures$interested_features_dataframe, file, row.names = FALSE)}
    )
    
    
    observeEvent(input$simulate_missingness, {
        
        show_modal_spinner(text = "Selecting Data")
        
        missingness <- summary(interaction(lapply(alignedFeatures$interested_features_dataframe[ , grepl( "mz" , names(alignedFeatures$interested_features_dataframe))],is.na)))
        
        missingness <- missingness[missingness>0]
        
        missingness_identified <- t(as.data.frame(strsplit(names(missingness), "[.]")))
        colnames(missingness_identified) <- as.matrix(strsplit(colnames(alignedFeatures$interested_features_dataframe[ , grepl( "mz" , names(alignedFeatures$interested_features_dataframe))]), "mz"))
        missingness_identified <- ifelse(missingness_identified=="TRUE", 1, 0)
        missingness_identified <- cbind(missingness_identified, num_missing = rowSums(missingness_identified))
        missingness_identified <- cbind(missingness_identified, missingness)
        missingness_identified <- cbind(missingness_identified, simulated_missingness = round(as.data.frame(missingness_identified)$missingness/sum(as.data.frame(missingness_identified)$missingness)*as.data.frame(missingness_identified)$missingness[which(as.data.frame(missingness_identified)$num_missing==0)],0))
        
        missingness_identified <- as.data.frame(missingness_identified)
        
        if(sum(missingness_identified$simulated_missingness)>missingness_identified$missingness[1])
        {
            missingness_identified$simulated_missingness[1] <- missingness_identified$simulated_missingness[1]-(sum(missingness_identified$simulated_missingness)-missingness_identified$missingness[1])
        }
        
        
        
        no_missingdata <<- alignedFeatures$interested_features_dataframe[rowSums(is.na(alignedFeatures$interested_features_dataframe[ , grepl( "mz" , names(alignedFeatures$interested_features_dataframe))])) == 0,] 
        
        alignedFeatures$nomissingdata <- no_missingdata
        
        no_missingdata[no_missingdata == 0] <- NA
        
        imputation_methods <- c("sv", "mn", "md", "bpca", "knn", "rf")
        
        no_missingdata_modeled_imputed <- NULL
        pca_variance <- list()
        
        remove_modal_spinner()
        
        nrmse <- NULL
        
        show_modal_progress_line()
        
        for(m in 1:length(imputation_methods))
        {
             no_missingdata_imputed <- mv_imputation(no_missingdata[,-1:-(length(complete_data$data)*2)], method = imputation_methods[m], check_df = FALSE)
            #no_missingdata_imputed <- mv_imputation(no_missingdata[,-1:-(length(test_modified)*2)], method = imputation_methods[m], check_df = FALSE)
            no_missingdata_imputed_copy <- no_missingdata_imputed
            
            no_missingdata_modified <- NULL
            
            for(i in 1:nrow(missingness_identified))
            {
                if(missingness_identified$simulated_missingness[i]>0)
                {
                    if(missingness_identified$num_missing[i]==0)
                    {
                        rand_ind <- sample(nrow(no_missingdata_imputed_copy), missingness_identified$simulated_missingness[i], replace = FALSE)
                        no_missingdata_modified <- rbind(no_missingdata_modified, no_missingdata_imputed_copy[rand_ind,])
                        no_missingdata_imputed_copy <- no_missingdata_imputed_copy[-rand_ind,]
                    }
                    else
                    {
                        rand_ind <- sample(nrow(no_missingdata_imputed_copy), missingness_identified$simulated_missingness[i], replace = FALSE)
                        selected_data <- no_missingdata_imputed_copy[rand_ind,]
                        missing_batches <- colnames(missingness_identified[which(missingness_identified[i,-(ncol(missingness_identified)-2):-ncol(missingness_identified)]==1)])
                        colnames_toedit <- foreach(j=1:length(missing_batches), .combine='c') %do% names(no_missingdata_imputed_copy[ , grepl( paste0(colnames(missingness_identified[which(missingness_identified[i,]==1)])[j], "[.]") , names( no_missingdata_imputed_copy ) ) ])
                        
                        for (k in 1:ncol(selected_data))
                        {
                            if (any(colnames(selected_data)[k]==colnames_toedit)==TRUE)
                            {
                                selected_data[,k] <- NA
                            }
                        }
                        
                        no_missingdata_modified <- rbind(no_missingdata_modified, selected_data)
                        no_missingdata_imputed_copy <- no_missingdata_imputed_copy[-rand_ind,] 
                    }
                }
                else next
            }
            
            
            no_missingdata_modified_imputed <- mv_imputation(no_missingdata_modified, method = imputation_methods[m], check_df = FALSE)
            
            
            no_missingdata_imputed_transformed <- as.matrix(log(no_missingdata_imputed, 2))
            no_missingdata_modified_imputed_transformed <- as.matrix(log(no_missingdata_modified_imputed, 2))
            
            na_vals <- which(is.na(no_missingdata_modified), arr.ind=TRUE)
            
            true <- NULL
            imp <- NULL
            
            for (i in 1:nrow(na_vals))
            {
                true <- c(true, no_missingdata_imputed_transformed[na_vals[i,1],na_vals[i,2]])
                imp <- c(imp, no_missingdata_modified_imputed_transformed[na_vals[i,1],na_vals[i,2]])
            }
            
            nrmse[m] <- rmse(true, imp, na.rm = FALSE)/mean(true)
            
            prin_comp <- prcomp(t(no_missingdata_modified_imputed_transformed), rank. = 2)
            components <- prin_comp[["x"]]
            components <- data.frame(components)
            explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
            explained_variance_ratio <- 100 * explained_variance_ratio
            components <- cbind(components, Batch = unlist(str_split(colnames(no_missingdata_imputed_transformed), "[.]", n = 2))[seq(1, length(unlist(str_split(colnames(no_missingdata_imputed_transformed), "[.]", n = 2))), 2)])
            components$PC2 <- -components$PC2
            components <- cbind(components, Method = imputation_methods[m])
            
            
            no_missingdata_modeled_imputed <- rbind(no_missingdata_modeled_imputed, components)
            pca_variance[[m]] <- explained_variance_ratio
            
            update_modal_progress(value = m/length(imputation_methods), text = paste("Peforming",imputation_methods[m], "imputation"))
            
        }
        remove_modal_progress()
        
        show_modal_spinner(text = "Generating Plots")
        output$modeledMissingnessplots <- renderPlotly({
            
            no_missingdata_modeled_imputed$Method <- factor(no_missingdata_modeled_imputed$Method, levels = c("sv", "mn",  "md", "bpca", "knn", "rf"))
            
            no_missingdata_modeled_imputed$Batch <- factor(no_missingdata_modeled_imputed$Batch, levels = names(complete_data$data))
            
            sv <- plot_ly(no_missingdata_modeled_imputed[no_missingdata_modeled_imputed$Method=="sv",], x = ~PC1, y = ~PC2, color = ~Batch, colors = hue_pal()(length(unique(no_missingdata_modeled_imputed$Batch))), type = 'scatter', mode = 'markers', legendgroup=~Batch) %>%
                layout(xaxis = list(title = paste('PC 1 (', round(pca_variance[[1]][1], 2), '%)', sep = '')), 
                       yaxis = list(title = paste('PC 1 (', round(pca_variance[[1]][2], 2), '%)'), sep = '')
                       )
            
            mn <- plot_ly(no_missingdata_modeled_imputed[no_missingdata_modeled_imputed$Method=="mn",], x = ~PC1, y = ~PC2, color = ~Batch, colors = hue_pal()(length(unique(no_missingdata_modeled_imputed$Batch))), type = 'scatter', mode = 'markers', legendgroup=~Batch) %>%
                layout(xaxis = list(title = paste('PC 1 (', round(pca_variance[[2]][1], 2), '%)', sep = '')), 
                       yaxis = list(title = paste('PC 1 (', round(pca_variance[[2]][2], 2), '%)'), sep = '')
                )
            
            md <- plot_ly(no_missingdata_modeled_imputed[no_missingdata_modeled_imputed$Method=="md",], x = ~PC1, y = ~PC2, color = ~Batch, colors = hue_pal()(length(unique(no_missingdata_modeled_imputed$Batch))), type = 'scatter', mode = 'markers', legendgroup=~Batch) %>%
                layout(xaxis = list(title = paste('PC 1 (', round(pca_variance[[3]][1], 2), '%)', sep = '')), 
                       yaxis = list(title = paste('PC 1 (', round(pca_variance[[3]][2], 2), '%)'), sep = '')
                )
            
            bpca <- plot_ly(no_missingdata_modeled_imputed[no_missingdata_modeled_imputed$Method=="bpca",], x = ~PC1, y = ~PC2, color = ~Batch, colors = hue_pal()(length(unique(no_missingdata_modeled_imputed$Batch))), type = 'scatter', mode = 'markers', legendgroup=~Batch) %>%
                layout(xaxis = list(title = paste('PC 1 (', round(pca_variance[[4]][1], 2), '%)', sep = '')), 
                       yaxis = list(title = paste('PC 1 (', round(pca_variance[[4]][2], 2), '%)'), sep = '')
                )
            
            knn <- plot_ly(no_missingdata_modeled_imputed[no_missingdata_modeled_imputed$Method=="knn",], x = ~PC1, y = ~PC2, color = ~Batch, colors = hue_pal()(length(unique(no_missingdata_modeled_imputed$Batch))), type = 'scatter', mode = 'markers', legendgroup=~Batch) %>%
                layout(xaxis = list(title = paste('PC 1 (', round(pca_variance[[5]][1], 2), '%)', sep = '')), 
                       yaxis = list(title = paste('PC 1 (', round(pca_variance[[5]][2], 2), '%)'), sep = '')
                )
            
            rf <- plot_ly(no_missingdata_modeled_imputed[no_missingdata_modeled_imputed$Method=="rf",], x = ~PC1, y = ~PC2, color = ~Batch, colors = hue_pal()(length(unique(no_missingdata_modeled_imputed$Batch))), type = 'scatter', mode = 'markers', legendgroup=~Batch) %>%
                layout(xaxis = list(title = paste('PC 1 (', round(pca_variance[[6]][1], 2), '%)', sep = '')), 
                       yaxis = list(title = paste('PC 1 (', round(pca_variance[[6]][2], 2), '%)'), sep = '')
                )
            
            subplot(style(sv, showlegend = F), style(mn, showlegend = F), style(md, showlegend = F), style(bpca, showlegend = F), style(knn, showlegend = F), rf, nrows = 2, titleY = TRUE, titleX = TRUE, margin = 0.07) %>%
            layout(plot_bgcolor='#e5ecf6',
                       annotations = list(
                            list(x = 0 , y = 1, text = paste( "Small Value; NRMSE = ", round(nrmse[1], 2)), showarrow = F, xref = "paper", yref= "paper", yanchor = "bottom"),
                            list(x = 0.5 , y = 1, text = paste( "Mean; NRMSE = ", round(nrmse[2], 2)), showarrow = F, xref = "paper", yref= "paper", yanchor = "bottom"),
                            list(x = 1 , y = 1, text = paste( "Median; NRMSE = ", round(nrmse[3], 2)), showarrow = F, xref = "paper", yref= "paper", yanchor = "bottom"),
                            list(x = 0, y = 0.42, text = paste("Bayesian PCA; NRMSE = ", round(nrmse[4], 2)), showarrow = F, xref = "paper", yref= "paper", yanchor = "bottom"),
                            list(x = 0.5, y = 0.42,text = paste("kNN; NRMSE = ", round(nrmse[5], 2)), showarrow = F, xref = "paper", yref= "paper", yanchor = "bottom"),
                            list(x = 1 , y = 0.42, text = paste("Random Forest; NRMSE = ", round(nrmse[6], 2)), showarrow = F, xref = "paper", yref= "paper", yanchor = "bottom")
                            ),
                       legend=list(title=list(text='Batch'))
            )

            })
        remove_modal_spinner()
        
    }, ignoreInit = TRUE)
    
    imputationMethod <- reactive({
        switch(input$impute,
               "small value" = "sv", 
               "mean" = "mn", 
               "median" = "md", 
               "bayesian PCA" = "bpca", 
               "kNN" = "knn", 
               "Random Forest" = "rf")
    })
    
    workingdatasets <- reactiveValues( imputed_dataset = NULL, waveica_batcheffect_corrected_dataset = NULL, combat_batcheffect_corrected_dataset = NULL)
    
    observeEvent(input$impute_data, {
        
        show_modal_spinner(text = "Imputing Data")
        
        output$NoMissingTrue <- renderPlotly({
            
            nomissing_data_PCA <- alignedFeatures$nomissingdata[,-1:-(length(complete_data$data)*2)]
            nomissing_data_PCA[nomissing_data_PCA == 0] <- NA
            nomissing_data_PCA <- mv_imputation(nomissing_data_PCA, method = imputationMethod() , check_df = FALSE)
            nomissing_data_PCA <- as.matrix(log(nomissing_data_PCA, 2))
            
            prin_comp_true <- prcomp(t(nomissing_data_PCA), rank. = 2)
            components_true <- prin_comp_true[["x"]]
            components_true <- data.frame(components_true)
            explained_variance_ratio <- summary(prin_comp_true)[["importance"]]['Proportion of Variance',]
            explained_variance_ratio <- 100 * explained_variance_ratio
            components_true <- cbind(components_true, Batch = unlist(str_split(colnames(nomissing_data_PCA), "[.]", n = 2))[seq(1, length(unlist(str_split(colnames(nomissing_data_PCA), "[.]", n = 2))), 2)])
            components_true$PC2 <- -components_true$PC2
            
            components_true$Batch <- factor(components_true$Batch, levels = names(complete_data$data))
            
            plot_ly(components_true, x = ~PC1, y = ~PC2, color = ~Batch, colors = hue_pal()(length(unique(components_true$Batch))), type = 'scatter', mode = 'markers', legendgroup=~Batch) %>%
                layout(
                    plot_bgcolor='#e5ecf6',
                    title = "Ground Truth of 'No Missing' Subset",
                    xaxis = list(title = paste('PC 1 (',toString(round(explained_variance_ratio[1],1)),'%)',sep = '')), 
                    yaxis = list(title = paste('PC 2 (',toString(round(explained_variance_ratio[2],1)),'%)',sep = '')),
                    legend=list(title=list(text='Batch'))
                )
            
        })
        
        entiredataset <- alignedFeatures$interested_features_dataframe
        
        entiredataset[entiredataset == 0] <- NA
        
        entiredataset_imputed <- mv_imputation(entiredataset[,-1:-(length(complete_data$data)*2)], method = imputationMethod(), check_df = FALSE)
        
        feature_info <- alignedFeatures$interested_features_dataframe[,1:(length(complete_data$data)*2)][ , match( c(paste0(refbatch_variables$refBatch_value,"mz"), paste0(refbatch_variables$refBatch_value,"rt")), names(alignedFeatures$interested_features_dataframe[1:(length(complete_data$data)*2)]))]
        
        entiredataset_imputed <- cbind(feature_info, entiredataset_imputed)
        
        entiredataset_imputed <<- entiredataset_imputed
        
        workingdatasets$imputed_dataset <- entiredataset_imputed
        
        output$ImputedData_Plot <- output$ImputedData_Plot2 <- renderPlotly({
            
            entiredataset_imputed_PCA <- as.matrix(log(entiredataset_imputed[,-1:-2], 2))
            
            prin_comp_entire <- prcomp(t(entiredataset_imputed_PCA), rank. = 2)
            components_entire <- prin_comp_entire[["x"]]
            components_entire <- data.frame(components_entire)
            explained_variance_ratio <- summary(prin_comp_entire)[["importance"]]['Proportion of Variance',]
            explained_variance_ratio <- 100 * explained_variance_ratio
            components_entire <- cbind(components_entire, Batch = unlist(str_split(colnames(entiredataset_imputed_PCA), "[.]", n = 2))[seq(1, length(unlist(str_split(colnames(entiredataset_imputed_PCA), "[.]", n = 2))), 2)])
            components_entire$PC2 <- -components_entire$PC2
            
            components_entire$Batch <- factor(components_entire$Batch, levels = names(complete_data$data))
            
            plot_ly(components_entire, x = ~PC1, y = ~PC2, color = ~Batch, colors = hue_pal()(length(unique(components_entire$Batch))), type = 'scatter', mode = 'markers', legendgroup=~Batch) %>%
                layout(
                    plot_bgcolor='#e5ecf6',
                    title = "Enitre Dataset Imputed",
                    xaxis = list(title = paste('PC 1 (',toString(round(explained_variance_ratio[1],1)),'%)',sep = '')), 
                    yaxis = list(title = paste('PC 2 (',toString(round(explained_variance_ratio[2],1)),'%)',sep = '')),
                    legend=list(title=list(text='Batch'))
                )
            
        })
        
       remove_modal_spinner() 
        
    }, ignoreInit = TRUE)
    
    output$download_imputed_data <- downloadHandler(
        filename = function() {paste("imputed_data-", Sys.Date(), ".csv", sep = "")},
        content = function(file) {
            showModal(modalDialog("Dowloading Imputed Data", footer = NULL))
            on.exit(removeModal())
            write.csv(workingdatasets$imputed_dataset, file, row.names = FALSE)}
    )
    
    additional_info <- reactive({
        read.csv(input$additionalinfo_file$datapath, header = TRUE, sep = ",")
        #read.csv("C:/Users/meh9jb/Desktop/Gates/MultiBatchPipeline/Gradient.csv", header = TRUE, sep = ",")
    })
    
    observeEvent(input$CorrectBatchEffect, {
        
        show_modal_spinner(text = "Performing WaveICA Correction")
        
        dataset <- t(workingdatasets$imputed_dataset[,-1:-2])
        
        colnames(dataset) <- paste(round(workingdatasets$imputed_dataset[,2],2), round(workingdatasets$imputed_dataset[,1],4) ,sep = "_")
        
        dataset <- dataset[order(match(sapply(strsplit(rownames(dataset), "[.]"), "[", 2), additional_info()[,1])),]
        
        dataset_WaveICA_corrected <- WaveICA_2.0(dataset, Injection_Order = additional_info()[,4], alpha = 0, Cutoff = 0.1, K=10)
        
        workingdatasets$waveica_batcheffect_corrected_dataset <- t(dataset_WaveICA_corrected$data_wave)
        
        output$WaveICACorrected <- renderPlotly({
            
            dataset_PCA <- scale(dataset_WaveICA_corrected$data_wave)
            
            prin_comp_entire <- prcomp(dataset_PCA, rank. = 2)
            components_entire <- prin_comp_entire[["x"]]
            components_entire <- data.frame(components_entire)
            explained_variance_ratio <- summary(prin_comp_entire)[["importance"]]['Proportion of Variance',]
            explained_variance_ratio <- 100 * explained_variance_ratio
            components_entire <- cbind(components_entire, Batch = unlist(str_split(rownames(dataset_WaveICA_corrected$data_wave), "[.]", n = 2))[seq(1, length(unlist(str_split(rownames(dataset_WaveICA_corrected$data_wave), "[.]", n = 2))), 2)])
            components_entire <- cbind(components_entire, SampleType = additional_info()[,3][match(sapply(strsplit(rownames(components_entire), "[.]"), "[", 2), additional_info()[,1])])
            components_entire <- cbind(components_entire, SampleName = additional_info()[,1][match(sapply(strsplit(rownames(components_entire), "[.]"), "[", 2), additional_info()[,1])])
            components_entire$PC2 <- -components_entire$PC2
            
            components_entire$Batch <- factor(components_entire$Batch, levels = names(complete_data$data))
            
            components_entire %>% 
                plot_ly(x=~PC1,
                        y=~PC2,
                        type="scatter", 
                        mode="markers",
                        colors = hue_pal()(length(unique(components_entire$Batch))), 
                        color=~Batch, 
                        symbol=~SampleType, 
                        legendgroup=~SampleType,
                        text = ~SampleName,
                        marker = list(size = 10)) %>%
                layout(
                    title = "Batch Effect Corrected: WaveICA 2.0",
                    xaxis = list(title = paste('PC 1 (',toString(round(explained_variance_ratio[1],1)),'%)',sep = '')), 
                    yaxis = list(title = paste('PC 2 (',toString(round(explained_variance_ratio[2],1)),'%)',sep = '')),
                    legend=list(title=list(text='Batch'))
                )
            
        })
        
        update_modal_spinner(text = "Performing ComBat Correction")
        
        dataset_combat_corrected <- ComBat(t(dataset), batch = sapply(strsplit(rownames(dataset), "[.]"), "[", 1), mod=NULL)
      
        workingdatasets$combat_batcheffect_corrected_dataset <- dataset_combat_corrected
        
        output$ComBatCorrected <- renderPlotly({
            
            dataset_PCA <- scale(t(dataset_combat_corrected))
            
            prin_comp_entire <- prcomp(dataset_PCA, rank. = 2)
            components_entire <- prin_comp_entire[["x"]]
            components_entire <- data.frame(components_entire)
            explained_variance_ratio <- summary(prin_comp_entire)[["importance"]]['Proportion of Variance',]
            explained_variance_ratio <- 100 * explained_variance_ratio
            components_entire <- cbind(components_entire, Batch = unlist(str_split(colnames(dataset_combat_corrected), "[.]", n = 2))[seq(1, length(unlist(str_split(colnames(dataset_combat_corrected), "[.]", n = 2))), 2)])
            components_entire <- cbind(components_entire, SampleType = additional_info()[,3][match(sapply(strsplit(rownames(components_entire), "[.]"), "[", 2), additional_info()[,1])])
            components_entire <- cbind(components_entire, SampleName = additional_info()[,1][match(sapply(strsplit(rownames(components_entire), "[.]"), "[", 2), additional_info()[,1])])
            components_entire$PC2 <- -components_entire$PC2
            
            components_entire$Batch <- factor(components_entire$Batch, levels = names(complete_data$data))
            
            components_entire %>% 
                plot_ly(x=~PC1,
                        y=~PC2,
                        type="scatter", 
                        mode="markers",
                        colors = hue_pal()(length(unique(components_entire$Batch))), 
                        color=~Batch, 
                        symbol=~SampleType, 
                        legendgroup=~SampleType,
                        text = ~SampleName,
                        marker = list(size = 10)) %>%
                layout(
                    title = "Batch Effect Corrected: ComBat",
                    xaxis = list(title = paste('PC 1 (',toString(round(explained_variance_ratio[1],1)),'%)',sep = '')), 
                    yaxis = list(title = paste('PC 2 (',toString(round(explained_variance_ratio[2],1)),'%)',sep = '')),
                    legend=list(title=list(text='Batch'))
                )
            
        })
        
        remove_modal_spinner()
        
    }, ignoreInit = TRUE)
    
    output$download_waveica_data <- downloadHandler(
        filename = function() {paste("waveica_corrected_data-", Sys.Date(), ".csv", sep = "")},
        content = function(file) {
            showModal(modalDialog("Dowloading WaveICA Corrected Data", footer = NULL))
            on.exit(removeModal())
            write.csv(workingdatasets$waveica_batcheffect_corrected_dataset, file)}
    )
    
    output$download_combat_data <- downloadHandler(
        filename = function() {paste("combat_corrected_data-", Sys.Date(), ".csv", sep = "")},
        content = function(file) {
            showModal(modalDialog("Dowloading ComBat Corrected Data", footer = NULL))
            on.exit(removeModal())
            write.csv(workingdatasets$combat_batcheffect_corrected_dataset, file)}
    )

}

# Run the application 
shinyApp(ui = ui, server = server)
