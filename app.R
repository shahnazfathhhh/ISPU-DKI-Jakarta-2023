# Import Library
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(readxl)
library(corrplot)
library(bslib)
library(shinycssloaders)
library(scales)
library(lubridate)
library(car)
library(zoo)
library(naniar)
library(stringr)
library(moments)
library(lmtest)     
library(AER)        
library(MASS)     
library(psych)      
library(dplyr)
library(tidyr)

#========
# DATA 
#========

df_2023 <- read_excel("data/fpevd.xlsx")

#======================
# DATA PREPROCESSING 
#======================

# Convert date to Date format
df_2023$tanggal <- as.Date(df_2023$tanggal)

# Convert "-" and "---" to NA
df_2023$pm_sepuluh[df_2023$pm_sepuluh %in% c("-", "---")] <- NA
df_2023$pm_sepuluh <- as.numeric(df_2023$pm_sepuluh)
df_2023$pm_duakomalima[df_2023$pm_duakomalima %in% c("-", "---")] <- NA
df_2023$pm_duakomalima <- as.numeric(df_2023$pm_duakomalima)
df_2023$sulfur_dioksida[df_2023$sulfur_dioksida %in% c("-", "---")] <- NA
df_2023$sulfur_dioksida <- as.numeric(df_2023$sulfur_dioksida)
df_2023$karbon_monoksida[df_2023$karbon_monoksida %in% c("-", "---")] <- NA
df_2023$karbon_monoksida <- as.numeric(df_2023$karbon_monoksida)
df_2023$ozon[df_2023$ozon %in% c("-", "---")] <- NA
df_2023$ozon <- as.numeric(df_2023$ozon)
df_2023$nitrogen_dioksida[df_2023$nitrogen_dioksida %in% c("-", "---")] <- NA
df_2023$nitrogen_dioksida <- as.numeric(df_2023$nitrogen_dioksida)
df_2023$max[df_2023$max %in% c("-", "---")] <- NA
df_2023$max <- as.numeric(df_2023$max)

df_2023$stasiun <- as.factor(df_2023$stasiun)

# raw_data SEBELUM preprocessing
raw_data <- df_2023  

# Missing values summary BEFORE cleaning
parameter <- c("pm_sepuluh", "pm_duakomalima", "sulfur_dioksida",
               "karbon_monoksida", "ozon", "nitrogen_dioksida")

missing_before <- data.frame(
  Variable = names(raw_data),
  Missing_Count = colSums(is.na(raw_data)),
  Missing_Percent = round(colSums(is.na(raw_data)) / nrow(raw_data) * 100, 2)
)

# ===============
# PREPROCESSING 
# ===============

df_cleaned <- df_2023

# 1. Interpolation 
df_cleaned <- df_cleaned %>%
  arrange(tanggal, stasiun) %>%
  group_by(stasiun) %>%
  mutate(
    pm_sepuluh = round(zoo::na.approx(pm_sepuluh, na.rm = FALSE)),
    pm_duakomalima = round(zoo::na.approx(pm_duakomalima, na.rm = FALSE)),
    sulfur_dioksida = round(zoo::na.approx(sulfur_dioksida, na.rm = FALSE)),
    karbon_monoksida = round(zoo::na.approx(karbon_monoksida, na.rm = FALSE)),
    ozon = round(zoo::na.approx(ozon, na.rm = FALSE)),
    nitrogen_dioksida = round(zoo::na.approx(nitrogen_dioksida, na.rm = FALSE))
  ) %>%
  ungroup()

# 2. Median imputation 
df_cleaned <- df_cleaned %>%
  group_by(tanggal) %>%
  mutate(across(all_of(parameter), ~ ifelse(is.na(.x),
                                            round(median(.x, na.rm = TRUE)),
                                            .x))) %>%
  ungroup()

# 3. Calculate max ISPU 
df_cleaned <- df_cleaned %>%
  rowwise() %>%
  mutate(
    max = max(c_across(all_of(parameter)), na.rm = TRUE),
    parameter_pencemar_kritis = parameter[which.max(c_across(all_of(parameter)))]
  ) %>%
  ungroup()

# 4. Categorize ISPU
df_cleaned <- df_cleaned %>%
  mutate(
    kategori = case_when(
      max >= 0 & max <= 50 ~ "BAIK",
      max >= 51 & max <= 100 ~ "SEDANG",
      max >= 101 & max <= 199 ~ "TIDAK SEHAT",
      max >= 200 & max <= 299 ~ "SANGAT TIDAK SEHAT",
      max >= 300 ~ "BERBAHAYA",
      TRUE ~ "TIDAK ADA DATA"
    )
  )

# 5. Convert to factors
df_cleaned$kategori <- as.factor(df_cleaned$kategori)
df_cleaned$parameter_pencemar_kritis <- as.factor(df_cleaned$parameter_pencemar_kritis)

# 6. Handle name inconsistencies
df_cleaned <- df_cleaned %>%
  mutate(
    stasiun = str_squish(stasiun),
    stasiun = trimws(stasiun),
    stasiun = case_when(
      stasiun %in% c("DKI5 Kebon Jeruk Jakarta Barat", "DKI5 Kebon Jeruk") ~ "DKI5 Kebon Jeruk",
      TRUE ~ stasiun
    ),
    stasiun = as.factor(stasiun)
  )

# 7. Add temporal variables for regression
df_cleaned$tahun <- as.numeric(format(df_cleaned$tanggal, "%Y"))
df_cleaned$bulan <- as.numeric(format(df_cleaned$tanggal, "%m"))
df_cleaned$hari <- as.numeric(format(df_cleaned$tanggal, "%d"))
df_cleaned$bulan_nama <- format(df_cleaned$tanggal, "%B")

# 8. Remove redundant column
df_cleaned <- df_cleaned %>%
  dplyr::select(-periode_data)

# Missing values summary AFTER cleaning
missing_after <- data.frame(
  Variable = names(df_cleaned),
  Missing_Count = colSums(is.na(df_cleaned)),
  Missing_Percent = round(colSums(is.na(df_cleaned)) / nrow(df_cleaned) * 100, 2)
)

# =============
# OUTLIER 
# =============

# Initial model
model1 <- lm(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida +
               karbon_monoksida + nitrogen_dioksida +
               stasiun + kategori + tahun + bulan + hari,
             data = df_cleaned)

# Calculate influential observations
hat_values <- hatvalues(model1)
student_resid <- rstudent(model1)
cooks_d <- cooks.distance(model1)

n <- nrow(df_cleaned)
p <- length(coef(model1)) - 1
h_bar <- (p + 1) / n
threshold_3h <- 3 * h_bar

# Identify observations to remove
high_leverage_index <- which(hat_values > threshold_3h)
high_student_resid_index <- which(abs(student_resid) > 2)
high_cooks_distance_index <- which(cooks_d > (4 / n))

remove_index <- unique(c(high_leverage_index,
                         high_student_resid_index,
                         high_cooks_distance_index))

# Create cleaned dataset without influential observations
df_model_clean <- df_cleaned[-remove_index, ]

# Build model without influential observations
model_clean <- lm(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida +
                    karbon_monoksida + nitrogen_dioksida +
                    stasiun + kategori + tahun + bulan + hari,
                  data = df_model_clean)

# Standardization 
df_standardized <- df_model_clean %>%
  mutate(
    pm10_std = scale(pm_sepuluh)[,1],
    pm2.5_std = scale(pm_duakomalima)[,1],
    so2_std = scale(sulfur_dioksida)[,1],
    co_std = scale(karbon_monoksida)[,1],
    o3_std = scale(ozon)[,1],
    no2_std = scale(nitrogen_dioksida)[,1]
  )

model_standardized <- lm(o3_std ~ pm10_std + pm2.5_std + so2_std + co_std + no2_std +
                           stasiun + kategori + tahun + bulan + hari,
                         data = df_standardized)

# Use cleaned data for analysis 
analysis_data <- df_model_clean

# ============================================================================
# UI DEFINITION
# ============================================================================
ui <- dashboardPage(
  skin = "green",
  
  dashboardHeader(
    title = "ISPU DKI Jakarta 2023",
    titleWidth = 300
  ),
  
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("Introduction", tabName = "intro", icon = icon("home")),
      menuItem("Data Sources", tabName = "data_sources", icon = icon("database")),
      menuItem("Data Preprocessing", tabName = "preprocessing", icon = icon("broom")),
      menuItem("Exploratory Analysis", tabName = "eda", icon = icon("chart-bar")),
      menuItem("Regression Analysis", tabName = "regression", icon = icon("line-chart"))
    ),
    hr(),
    selectInput("filter_station", "Filter by Station:", 
                choices = c("All", levels(analysis_data$stasiun)),
                selected = "All"),
    selectInput("filter_month", "Filter by Month:", 
                choices = c("All", unique(analysis_data$bulan_nama)),
                selected = "All"),
    actionButton("reset_filters", "Reset Filters", 
                 icon = icon("redo"),
                 style = "margin: 10px 15px; width: 90%;")
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .box { border-top: 3px solid #00a65a; }
        .content-wrapper { background-color: #f4f6f9; }
        .sdg-box { 
          background: linear-gradient(135deg, #1e7e34 0%, #28a745 100%);
          color: white;
          padding: 20px;
          border-radius: 10px;
          margin-bottom: 20px;
        }
        .info-card {
          background: white;
          padding: 15px;
          border-radius: 8px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
          margin: 10px 0;
        }
      "))
    ),
    
    tabItems(
      # ======================
      # TAB 1: INTRODUCTION
      # ======================
      tabItem(
        tabName = "intro",
        fluidRow(
          box(
            width = 12,
            solidHeader = TRUE,
            div(class = "sdg-box",
                h2(icon("leaf"), " Analisis Indeks Standar Pencemaran Udara (ISPU) DKI Jakarta 2023"),
                h3("Kelompok 9"),
                h4("Supporting SDG 13: Climate Action")
            )
          )
        ),
        fluidRow(
          box(
            title = "Latar Belakang",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            div(class = "info-card",
                p("Pencemaran udara menjadi ancaman serius bagi lingkungan dan kesehatan manusia, di mana emisi gas rumah kaca dan polutan seperti CO, NO₂, SO₂, PM₁₀, dan O₃ turut mempercepat perubahan iklim serta meningkatkan risiko penyakit pernapasan hingga kematian dini secara global. Jakarta sebagai pusat aktivitas ekonomi menghadapi beban polusi tinggi yang didominasi oleh sektor transportasi dan industri, dengan konsentrasi polutan yang sering melampaui standar aman WHO. Meskipun pemantauan kualitas udara melalui Indeks Standar Pencemar Udara (ISPU) telah dilakukan sejak 2010, analisis mendalam masih dibutuhkan untuk memahami hubungan antarpolutan, terutama peran prekursor dalam pembentukan ozon troposferik. Penelitian ini menerapkan regresi binomial negatif untuk mengidentifikasi variabel yang memengaruhi O₃, sehingga dapat mendukung kebijakan pengendalian emisi dan strategi mitigasi perubahan iklim di Jakarta."),
                h4("SDG 13: Climate Action"),
                p("Analisis ini mendukung Tujuan Pembangunan Berkelanjutan (SDG) ke-13 
                  yang berfokus pada pengambilan tindakan mendesak untuk memerangi perubahan 
                  iklim dan dampaknya. Memahami pola kualitas udara sangat penting untuk 
                  mengembangkan strategi aksi iklim yang efektif."),
                h4("Indeks Standar Pencemar Udara (ISPU)"),
                p("ISPU adalah angka yang digunakan untuk menggambarkan kondisi kualitas 
                  udara di suatu lokasi tertentu. Semakin tinggi nilai ISPU, semakin buruk 
                  kualitas udara dan semakin besar dampaknya terhadap kesehatan.")
            )
          ),
          box(
            title = "Rumusan Masalah",
            status = "warning",
            solidHeader = TRUE,
            width = 6,
            div(class = "info-card",
                h4("Pertanyaan Penelitian"),
                tags$ol(
                  tags$li("Bagaimana karakteristik pemusatan dan persebaran Indeks Standar Pencemaran Udara (ISPU) di berbagai stasiun pemantauan DKI Jakarta pada tahun 2023?"),
                  tags$li("Bagaimana model yang tepat untuk memodelkan data kualitas udara DKI Jakarta pada tahun 2023?"),
                  tags$li("Bagaimana faktor-faktor yang berpengaruh signifikan terhadap parameter Ozon di DKI Jakarta?")
                )
            )
          )
        ),
        fluidRow(
          box(
            title = "Tujuan Penelitian",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            div(class = "info-card",
                h4("Tujuan Utama"),
                tags$ol(
                  tags$li("Mengetahui karakteristik pemusatan dan persebaran Indeks Standar Pencemaran Udara (ISPU) di berbagai stasiun pemantauan DKI Jakarta pada tahun 2023."),
                  tags$li("Mengetahui model yang tepat untuk memodelkan data kualitas udara DKI Jakarta pada tahun 2023."),
                  tags$li("Mengetahui faktor-faktor yang berpengaruh signifikan terhadap parameter Ozon di DKI Jakarta."),
                )
            )
          )
        )
      ),
      
      # =========================
      # TAB 2: DATA SOURCES
      # =========================
      tabItem(
        tabName = "data_sources",
        h2("Sumber Data"),
        fluidRow(
          box(
            title = "Informasi Data",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            div(class = "info-card",
                h4(icon("database"), " Sumber"),
                p(strong("Portal Data:"), " Satu Data Jakarta"),
                p(strong("URL:"), 
                  a("https://satudata.jakarta.go.id", 
                    href = "https://satudata.jakarta.go.id/open-data/detail?kategori=dataset&page_url=data-indeks-standar-pencemar-udara-ispu-di-provinsi-dki-jakarta-2023&data_no=1",
                    target = "_blank")),
                p(strong("Dataset:"), "Indeks Standar Pencemar Udara (ISPU) - DKI Jakarta 2023"),
                p(strong("Format:"), "Excel (.xlsx)"),
                h4(icon("map-marker"), " Stasiun Pemantauan"),
                p("Data dikumpulkan dari 5 stasiun pemantauan di seluruh DKI Jakarta:",
                  tags$ul(
                    tags$li("DKI1 Bunderan HI (Jakarta Pusat)"),
                    tags$li("DKI2 Kelapa Gading (Jakarta Utara)"),
                    tags$li("DKI3 Jagakarsa (Jakarta Selatan)"),
                    tags$li("DKI4 Lubang Buaya (Jakarta Timur)"),
                    tags$li("DKI5 Kebon Jeruk (Jakarta Barat)")
                  ))
            )
          )
        ),
        fluidRow(
          box(
            title = "Deskripsi Variabel",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            DTOutput("variable_description")
          ),
          box(
            title = "Struktur Data",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            verbatimTextOutput("data_structure")
          )
        ),
        fluidRow(
          box(
            title = "Data Mentah",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            downloadButton("download_raw", "Download Data Mentah"),
            hr(),
            DTOutput("raw_data_table") %>% withSpinner()
          )
        )
      ),
      
      # ================================
      # TAB 3: DATA PREPROCESSING
      # ================================
      tabItem(
        tabName = "preprocessing",
        h2("Preprocessing Data"),
        
        fluidRow(
          box(
            title = "Tahapan Preprocessing",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            div(class = "info-card",
                h4("Langkah-langkah yang Diterapkan:"),
                tags$ol(
                  tags$li(strong("Konversi Tipe Data:"), " Mengubah tanggal ke format Date dan numerik ke numeric"),
                  tags$li(strong("Interpolasi Missing Values:"), " Mengisi nilai hilang menggunakan interpolasi linear per stasiun"),
                  tags$li(strong("Imputasi Median:"), " Mengisi nilai hilang yang tersisa dengan median per tanggal"),
                  tags$li(strong("Perhitungan ISPU:"), " Menghitung nilai maksimum dari semua polutan"),
                  tags$li(strong("Kategorisasi:"), " Mengklasifikasikan ISPU ke dalam kategori"),
                  tags$li(strong("Standarisasi Nama Stasiun:"), " Merapikan inkonsistensi nama stasiun"),
                  tags$li(strong("Deteksi & Penghapusan Influential Observations:"), " Menggunakan leverage, studentized residuals, dan Cook's distance"),
                  tags$li(strong("Standarisasi Data:"), " Z-score standardization untuk variabel numerik")
                )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Tipe Data",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            verbatimTextOutput("data_types")
          ),
          box(
            title = "Analisis Missing Values",
            status = "warning",
            solidHeader = TRUE,
            width = 6,
            h4("Sebelum Preprocessing:"),
            DTOutput("missing_before_table"),
            hr(),
            h4("Setelah Preprocessing:"),
            DTOutput("missing_after_table")
          )
        ),
        
        fluidRow(
          box(
            title = "Visualisasi Missing Values",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("missing_visualization", height = 500) %>% withSpinner()
          )
        ),
        
        fluidRow(
          box(
            title = "Outlier Box Plots",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            selectInput("outlier_var", "Pilih Variabel:",
                        choices = c("pm_sepuluh", "pm_duakomalima", "sulfur_dioksida",
                                    "karbon_monoksida", "ozon", "nitrogen_dioksida", "max")),
            plotlyOutput("outlier_boxplot", height = 400) %>% withSpinner()
          )
        ),
        
        fluidRow(
          box(
            title = "Influential Observations Analysis",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            uiOutput("influential_summary"),
            hr(),
            plotOutput("influence_plot", height = 600) %>% withSpinner()
          )
        ),
        
        fluidRow(
          box(
            title = "Perbandingan Model",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DTOutput("model_comparison_table")
          )
        ),
        
        fluidRow(
          box(
            title = "Data Bersih",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            downloadButton("download_clean", "Download Data Bersih"),
            hr(),
            DTOutput("clean_data_table") %>% withSpinner()
          )
        )
      ),
      
      # =======================================
      # TAB 4: EXPLORATORY DATA ANALYSIS
      # =======================================
      tabItem(
        tabName = "eda",
        h2("Analisis Data Eksploratif"),
        fluidRow(
          valueBoxOutput("avg_ispu", width = 3),
          valueBoxOutput("max_ispu", width = 3),
          valueBoxOutput("unhealthy_days", width = 3),
          valueBoxOutput("most_polluted_station", width = 3)
        ),
        fluidRow(
          box(
            title = "Ringkasan Statistik Komprehensif",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            DTOutput("comprehensive_stats")
          )
        ),
        fluidRow(
          box(
            title = "Distribusi ISPU Berdasarkan Kategori",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            plotlyOutput("ispu_category_pie", height = 350) %>% withSpinner()
          ),
          box(
            title = "Rata-rata ISPU per Stasiun",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            plotlyOutput("ispu_by_station", height = 350) %>% withSpinner()
          )
        ),
        fluidRow(
          box(
            title = "Tren ISPU Sepanjang Waktu",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("ispu_time_series", height = 400) %>% withSpinner()
          )
        ),
        fluidRow(
          box(
            title = "Perbandingan Konsentrasi Polutan",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("pollutant_comparison", height = 400) %>% withSpinner()
          )
        )
      ),
      
      # ================================
      # TAB 5: REGRESSION ANALYSIS
      # ================================
      tabItem(
        tabName = "regression",
        h2("Analisis Regresi: Prediksi Tingkat O₃"),
        
        fluidRow(
          box(
            title = "GLM Poisson",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            h4("Model GLM Poisson"),
            verbatimTextOutput("poisson_summary"),
            hr(),
            h4("Uji Dispersi"),
            verbatimTextOutput("dispersion_test"),
            uiOutput("dispersion_interpretation")
          )
        ),
        
        fluidRow(
          box(
            title = "Model GLM Binomial Negatif",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            h4("Model: O₃ ~ PM₁₀ + PM₂.₅ + SO₂ + CO + NO₂ + Stasiun + Kategori + Tahun + Bulan + Hari"),
            verbatimTextOutput("negbin_summary"),
            hr()
          )
        ),
        
        fluidRow(
          box(
            title = "Perbandingan Model",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            DTOutput("model_comparison"),
            hr(),
          ),
          box(
            title = "Performa Model Terpilih",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            uiOutput("best_model_performance")
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER LOGIC
# ============================================================================
server <- function(input, output, session) {
  
  # Reactive filtered data
  filtered_data <- reactive({
    data <- analysis_data
    
    if (input$filter_station != "All") {
      data <- data %>% filter(stasiun == input$filter_station)
    }
    
    if (input$filter_month != "All") {
      data <- data %>% filter(bulan_nama == input$filter_month)
    }
    
    data
  })
  
  # Reset filters
  observeEvent(input$reset_filters, {
    updateSelectInput(session, "filter_station", selected = "All")
    updateSelectInput(session, "filter_month", selected = "All")
  })
  
  # =================================
  # TAB 2: DATA SOURCES OUTPUTS
  # =================================
  output$variable_description <- renderDT({
    var_desc <- data.frame(
      Variabel = c("tanggal", "stasiun", "pm_sepuluh", "pm_duakomalima", 
                   "sulfur_dioksida", "karbon_monoksida", "ozon", 
                   "nitrogen_dioksida", "max", "parameter_pencemar_kritis", "kategori"),
      Deskripsi = c(
        "Tanggal pengukuran",
        "Lokasi stasiun pemantauan",
        "Particulate Matter 10 μm (μg/m³)",
        "Particulate Matter 2.5 μm (μg/m³)",
        "Sulfur Dioksida (μg/m³)",
        "Karbon Monoksida (μg/m³)",
        "Ozon troposfer (μg/m³)",
        "Nitrogen Dioksida (μg/m³)",
        "Nilai ISPU maksimum",
        "Parameter dengan nilai ISPU tertinggi",
        "Kategori kualitas udara"
      ),
      Tipe = c("Date", "Factor", rep("Numeric", 7), "Factor", "Factor")
    )
    datatable(var_desc, options = list(pageLength = 15, dom = 't'), rownames = FALSE)
  })
  
  output$data_structure <- renderPrint({
    str(raw_data)
  })
  
  output$raw_data_table <- renderDT({
    datatable(head(raw_data, 100), 
              options = list(pageLength = 10, scrollX = TRUE),
              filter = "top",
              rownames = FALSE)
  })
  
  output$download_raw <- downloadHandler(
    filename = function() {
      paste0("ISPU_Data_Mentah_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(raw_data, file, row.names = FALSE)
    }
  )
  
  # =================================
  # TAB 3: PREPROCESSING OUTPUTS
  # =================================
  
  output$data_types <- renderPrint({
    cat("Tipe Data Setelah Preprocessing:\n\n")
    sapply(df_cleaned, class)
  })
  
  output$missing_before_table <- renderDT({
    datatable(head(missing_before, 12), 
              options = list(pageLength = 12, dom = 't'),
              rownames = FALSE)
  })
  
  output$missing_after_table <- renderDT({
    datatable(head(missing_after, 12), 
              options = list(pageLength = 12, dom = 't'),
              rownames = FALSE)
  })
  
  output$missing_visualization <- renderPlot({
    par(mfrow = c(1, 2), mar = c(8, 4, 4, 2))
    
    missing_pct_before <- colSums(is.na(raw_data[, parameter])) / nrow(raw_data) * 100
    missing_pct_after <- colSums(is.na(df_cleaned[, parameter])) / nrow(df_cleaned) * 100
    
    max_y <- max(c(missing_pct_before, missing_pct_after), na.rm = TRUE) + 5
    
    barplot(missing_pct_before, 
            main = "Missing Values SEBELUM Preprocessing",
            ylab = "Persentase Missing (%)",
            las = 2,
            col = "#dd4b39",
            ylim = c(0, max_y),
            cex.names = 0.8)
    
    barplot(missing_pct_after, 
            main = "Missing Values SETELAH Preprocessing",
            ylab = "Persentase Missing (%)",
            las = 2,
            col = "#00a65a",
            ylim = c(0, max_y),
            cex.names = 0.8)
  })
  
  output$outlier_boxplot <- renderPlotly({
    plot_ly(df_cleaned, y = ~get(input$outlier_var), type = "box",
            name = input$outlier_var,
            marker = list(color = "#00a65a")) %>%
      layout(title = paste("Deteksi Outlier:", input$outlier_var),
             yaxis = list(title = input$outlier_var))
  })
  
  output$influential_summary <- renderUI({
    HTML(paste(
      "<div class='info-card'>",
      "<h4>Ringkasan Influential Observations</h4>",
      "<p><strong>Total observasi awal:</strong>", format(nrow(df_cleaned), big.mark = ","), "</p>",
      "<p><strong>High leverage (> 3h̄):</strong>", length(high_leverage_index), "</p>",
      "<p><strong>High studentized residuals (|t| > 2):</strong>", length(high_student_resid_index), "</p>",
      "<p><strong>High Cook's distance (> 4/n):</strong>", length(high_cooks_distance_index), "</p>",
      "<p><strong>Total observasi dihapus:</strong>", length(remove_index), 
      "(", round(length(remove_index)/nrow(df_cleaned)*100, 2), "%)</p>",
      "<p><strong>Observasi dalam model final:</strong>", format(nrow(analysis_data), big.mark = ","), "</p>",
      "</div>"
    ))
  })
  
  output$influence_plot <- renderPlot({
    par(mfrow = c(2, 2))
    
    plot(hat_values, main = "Leverage Values", 
         ylab = "Leverage", xlab = "Index",
         pch = 19, col = ifelse(hat_values > threshold_3h, "red", "blue"))
    abline(h = threshold_3h, col = "red", lty = 2, lwd = 2)
    
    plot(student_resid, main = "Studentized Residuals",
         ylab = "Studentized Residuals", xlab = "Index",
         pch = 19, col = ifelse(abs(student_resid) > 2, "red", "blue"))
    abline(h = c(-2, 2), col = "red", lty = 2, lwd = 2)
    
    plot(cooks_d, main = "Cook's Distance",
         ylab = "Cook's Distance", xlab = "Index",
         pch = 19, col = ifelse(cooks_d > (4/n), "red", "blue"))
    abline(h = 4/n, col = "red", lty = 2, lwd = 2)
    
    influencePlot(model1, main = "Influence Plot",
                  sub = "Ukuran lingkaran = Cook's distance")
  })
  
  output$model_comparison_table <- renderDT({
    comparison <- data.frame(
      Model = c(
        "Model tanpa penanganan outlier",
        "Model setelah menghapus influential observations",
        "Model setelah standardisasi"
      ),
      `Jumlah Observasi` = c(
        nrow(df_cleaned),
        nrow(df_model_clean),
        nrow(df_standardized)
      ),
      `R-squared` = c(
        summary(model1)$r.squared,
        summary(model_clean)$r.squared,
        summary(model_standardized)$r.squared
      ),
      `Adjusted R-squared` = c(
        summary(model1)$adj.r.squared,
        summary(model_clean)$adj.r.squared,
        summary(model_standardized)$adj.r.squared
      ),
      `RMSE` = c(
        sqrt(mean(model1$residuals^2)),
        sqrt(mean(model_clean$residuals^2)),
        sqrt(mean(model_standardized$residuals^2))
      ),
      check.names = FALSE
    )
    
    datatable(comparison, 
              options = list(pageLength = 5, dom = 't'),
              rownames = FALSE) %>%
      formatRound(columns = c("R-squared", "Adjusted R-squared", "RMSE"), digits = 4)
  })
  
  output$clean_data_table <- renderDT({
    datatable(head(analysis_data, 100), 
              options = list(pageLength = 10, scrollX = TRUE),
              filter = "top",
              rownames = FALSE)
  })
  
  output$download_clean <- downloadHandler(
    filename = function() {
      paste0("ISPU_Data_Bersih_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(analysis_data, file, row.names = FALSE)
    }
  )
  
  # ============================
  # TAB 4: EDA OUTPUTS
  # ============================
  
  output$avg_ispu <- renderValueBox({
    valueBox(
      round(mean(filtered_data()$max, na.rm = TRUE), 1),
      "Rata-rata ISPU",
      icon = icon("wind"),
      color = "aqua"
    )
  })
  
  output$max_ispu <- renderValueBox({
    valueBox(
      max(filtered_data()$max, na.rm = TRUE),
      "ISPU Maksimum",
      icon = icon("exclamation-triangle"),
      color = "red"
    )
  })
  
  output$unhealthy_days <- renderValueBox({
    unhealthy <- filtered_data() %>%
      filter(kategori %in% c("TIDAK SEHAT", "SANGAT TIDAK SEHAT", "BERBAHAYA")) %>%
      distinct(tanggal) %>%
      nrow()
    
    valueBox(
      unhealthy,
      "Hari Tidak Sehat",
      icon = icon("frown"),
      color = "orange"
    )
  })
  
  output$most_polluted_station <- renderValueBox({
    station <- filtered_data() %>%
      group_by(stasiun) %>%
      summarise(avg_ispu = mean(max, na.rm = TRUE)) %>%
      arrange(desc(avg_ispu)) %>%
      slice(1) %>%
      pull(stasiun) %>%
      as.character()
    
    valueBox(
      station,
      "Stasiun Terpolusi",
      icon = icon("map-pin"),
      color = "maroon"
    )
  })
  
  output$comprehensive_stats <- renderDT({
    stats <- filtered_data() %>%
      dplyr::select(pm_sepuluh, pm_duakomalima, sulfur_dioksida, 
                    karbon_monoksida, ozon, nitrogen_dioksida, max) %>%
      tidyr::pivot_longer(everything(), names_to = "Variable") %>%
      dplyr::mutate(Variable = dplyr::case_when(
        Variable == "pm_sepuluh" ~ "PM10",
        Variable == "pm_duakomalima" ~ "PM2.5",
        Variable == "sulfur_dioksida" ~ "SO2",
        Variable == "karbon_monoksida" ~ "CO",
        Variable == "ozon" ~ "O3",
        Variable == "nitrogen_dioksida" ~ "NO2",
        Variable == "max" ~ "ISPU",
        TRUE ~ Variable
      )) %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(
        Mean = round(mean(value, na.rm = TRUE), 2),
        Median = round(median(value, na.rm = TRUE), 2),
        SD = round(sd(value, na.rm = TRUE), 2),
        Min = round(min(value, na.rm = TRUE), 2),
        Max = round(max(value, na.rm = TRUE), 2),
        Q1 = round(quantile(value, 0.25, na.rm = TRUE), 2),
        Q3 = round(quantile(value, 0.75, na.rm = TRUE), 2)
      )
    
    DT::datatable(stats, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$ispu_category_pie <- renderPlotly({
    data <- filtered_data() %>%
      group_by(kategori) %>%
      summarise(Count = n()) %>%
      mutate(Percentage = round(Count/sum(Count)*100, 1))
    
    colors <- c("BAIK" = "#00a65a", "SEDANG" = "#3c8dbc", 
                "TIDAK SEHAT" = "#f39c12", "SANGAT TIDAK SEHAT" = "#dd4b39", 
                "BERBAHAYA" = "#000000")
    
    plot_ly(data, labels = ~kategori, values = ~Count, type = "pie",
            textinfo = "label+percent",
            marker = list(colors = colors[data$kategori])) %>%
      layout(title = "Distribusi Kategori ISPU")
  })
  
  output$ispu_by_station <- renderPlotly({
    data <- filtered_data() %>%
      group_by(stasiun) %>%
      summarise(Avg_ISPU = mean(max, na.rm = TRUE)) %>%
      arrange(desc(Avg_ISPU))
    
    plot_ly(data, x = ~stasiun, y = ~Avg_ISPU, type = "bar",
            marker = list(color = ~Avg_ISPU, colorscale = "Reds")) %>%
      layout(title = "Rata-rata ISPU per Stasiun",
             xaxis = list(title = "Stasiun"),
             yaxis = list(title = "Rata-rata ISPU"))
  })
  
  output$ispu_time_series <- renderPlotly({
    data <- filtered_data() %>%
      group_by(tanggal) %>%
      summarise(Avg_ISPU = mean(max, na.rm = TRUE))
    
    plot_ly(data, x = ~tanggal, y = ~Avg_ISPU, type = "scatter", mode = "lines",
            line = list(color = "#00a65a", width = 2)) %>%
      layout(title = "Tren ISPU Sepanjang Tahun 2023",
             xaxis = list(title = "Tanggal"),
             yaxis = list(title = "Rata-rata ISPU"))
  })
  
  output$pollutant_comparison <- renderPlotly({
    data <- filtered_data() %>%
      dplyr::select(pm_sepuluh, pm_duakomalima, sulfur_dioksida,
                    karbon_monoksida, ozon, nitrogen_dioksida) %>%
      tidyr::pivot_longer(everything(), 
                          names_to = "Pollutant", 
                          values_to = "Concentration") %>%
      dplyr::mutate(Pollutant = dplyr::case_when(
        Pollutant == "pm_sepuluh" ~ "PM10",
        Pollutant == "pm_duakomalima" ~ "PM2.5",
        Pollutant == "sulfur_dioksida" ~ "SO2",
        Pollutant == "karbon_monoksida" ~ "CO",
        Pollutant == "ozon" ~ "O3",
        Pollutant == "nitrogen_dioksida" ~ "NO2",
        TRUE ~ Pollutant
      ))
    
    plot_ly(data, x = ~Pollutant, y = ~Concentration, color = ~Pollutant,
            type = "box") %>%
      layout(title = "Distribusi Konsentrasi Polutan",
             yaxis = list(title = "Konsentrasi (μg/m³)"),
             showlegend = FALSE)
  })
  
  
  # ================================
  # TAB 5: REGRESSION OUTPUTS
  # ================================
  
  # GLM Poisson Model (SESUAI RMARKDOWN: tanpa lag)
  output$poisson_summary <- renderPrint({
    data <- filtered_data()
    model_poisson <- glm(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida +
                           karbon_monoksida + nitrogen_dioksida + stasiun + 
                           kategori + tahun + bulan + hari,
                         family = poisson(link = "log"),
                         data = data)
    summary(model_poisson)
  })
  
  output$dispersion_test <- renderPrint({
    data <- filtered_data()
    model_poisson <- glm(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida +
                           karbon_monoksida + nitrogen_dioksida + stasiun + 
                           kategori + tahun + bulan + hari,
                         family = poisson(link = "log"),
                         data = data)
    library(AER)
    dispersiontest(model_poisson, alternative = "greater")
  })
  
  output$dispersion_interpretation <- renderUI({
    data <- filtered_data()
    model_poisson <- glm(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida +
                           karbon_monoksida + nitrogen_dioksida + stasiun + 
                           kategori + tahun + bulan + hari,
                         family = poisson(link = "log"),
                         data = data)
    test <- AER::dispersiontest(model_poisson, alternative = "greater")
    
    if(test$p.value < 0.05) {
      HTML("<div class='alert alert-warning'><strong>⚠️ Overdispersion Detected:</strong> p-value < 0.05, 
           terjadi overdispersion. Model Binomial Negatif sesuai.</div>")
    } else {
      HTML("<div class='alert alert-success'><strong>✅ No Overdispersion:</strong> Model Poisson sesuai</div>")
    }
  })
  
  # Negative Binomial Model (SESUAI RMARKDOWN: tanpa lag)
  output$negbin_summary <- renderPrint({
    data <- filtered_data()
    
    library(MASS)
    model_negbin <- glm.nb(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida + 
                             karbon_monoksida + nitrogen_dioksida + stasiun + 
                             kategori + tahun + bulan + hari,
                           data = data)
    summary(model_negbin)
  })
  
  # Model Comparison
  output$model_comparison <- renderDT({
    data <- filtered_data()
    
    # Poisson Model
    model_poisson <- glm(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida +
                           karbon_monoksida + nitrogen_dioksida + stasiun + 
                           kategori + tahun + bulan + hari,
                         family = poisson(link = "log"),
                         data = data)
    
    # Negative Binomial Model
    library(MASS)
    model_negbin <- glm.nb(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida + 
                             karbon_monoksida + nitrogen_dioksida + stasiun + 
                             kategori + tahun + bulan + hari,
                           data = data)
    
    comparison <- data.frame(
      Model = c("GLM Poisson", "GLM Binomial Negatif"),
      AIC = c(AIC(model_poisson), AIC(model_negbin)),
      BIC = c(BIC(model_poisson), BIC(model_negbin)),
      Deviance = c(model_poisson$deviance, model_negbin$deviance)
    )
    
    datatable(comparison, 
              options = list(pageLength = 5, dom = 't'),
              rownames = FALSE) %>%
      formatRound(columns = c("AIC", "BIC", "Deviance"), digits = 2)
  })
  
  output$best_model_performance <- renderUI({
    data <- filtered_data()
    
    # Poisson Model
    model_poisson <- glm(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida +
                           karbon_monoksida + nitrogen_dioksida + stasiun + 
                           kategori + tahun + bulan + hari,
                         family = poisson(link = "log"),
                         data = data)
    
    # Negative Binomial Model
    library(MASS)
    model_negbin <- glm.nb(ozon ~ pm_sepuluh + pm_duakomalima + sulfur_dioksida + 
                             karbon_monoksida + nitrogen_dioksida + stasiun + 
                             kategori + tahun + bulan + hari,
                           data = data)
    
    aic_poisson <- AIC(model_poisson)
    aic_negbin <- AIC(model_negbin)
    
    best_model <- ifelse(aic_negbin < aic_poisson, "Binomial Negatif", "Poisson")
    best_aic <- min(aic_poisson, aic_negbin)
    
    HTML(paste(
      "<div class='info-card'>",
      "<h4>Model Terbaik</h4>",
      "<p><strong>Model Terpilih:</strong> GLM", best_model, "</p>",
      "<p><strong>AIC:</strong>", round(best_aic, 2), "</p>",
      "</div>"
    ))
  })
  
}

# ============================================================================
# RUN APPLICATION
# ============================================================================
shinyApp(ui = ui, server = server)