
ExtractDismissInformation_reshape <- function(path_to_pdf){
  library(pdftools)
  library(tidyverse)
  
  # load the pdf 
  testPDF <- pdf_text(path_to_pdf)
  
  ## subj name
  dat_name <- str_extract(testPDF, pattern = "姓名：.{2,3}")[1]
  dat_name2 <- str_split(dat_name, pattern = '：') %>% unlist() %>% .[2] %>% 
    trimws(which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  ## hospital ID
  dat_hospID <- str_extract(testPDF, pattern = "[0-9]{10}")[1]
  
  ## gender
  dat_gender <- str_extract(testPDF, pattern = "[男,女]")[1]
  dat_gender2 <- ifelse(dat_gender == '男', 'male', 'female')
  
  ## enroll time
  dat_enroll_date <- str_extract(testPDF, pattern = "入院时间：[0-9]{4}-[0-9]{2}-[0-9]{2}")[1]
  dat_enroll_date2 <- str_extract(dat_enroll_date, pattern = "[0-9]{4}-[0-9]{2}-[0-9]{2}")
  
  ## dismiss time
  dat_dismiss_date <- str_extract(testPDF, pattern = "出院时间：[0-9]{4}-[0-9]{2}-[0-9]{2}")[1]
  dat_dismiss_date2 <- str_extract(dat_dismiss_date, pattern = "[0-9]{4}-[0-9]{2}-[0-9]{2}")
  
  ## enroll diagnose
  dat_inhospital_diagnose <- str_extract(testPDF, pattern = "入院诊断：.*\n")[1]
  dat_enroll_diagnose <- ifelse(str_detect(dat_inhospital_diagnose, "重度抑郁"),"重度抑郁",
                                ifelse(str_detect(dat_inhospital_diagnose, "复发性抑郁"), "复发性抑郁",
                                       ifelse(str_detect(dat_inhospital_diagnose, "心境障碍"), "心境障碍",
                                              ifelse(str_detect(dat_inhospital_diagnose, "抑郁发作"), "抑郁发作",'其他'))))
  dat_enroll_psycho <- ifelse(str_detect(dat_inhospital_diagnose, "不伴精神病"),'No',
                              ifelse(str_detect(dat_inhospital_diagnose, "，伴精神病"), 'Yes',
                                     NA))
  dat_enroll_soma <- ifelse(str_detect(dat_inhospital_diagnose, "不伴躯体化"),'No',
                            ifelse(str_detect(dat_inhospital_diagnose, "，伴躯体化"), 'Yes',
                                   NA))
  ## dismiss diagnose
  dat_outhospital_diagnose <- str_extract(testPDF, pattern = "出院诊断：.*\n")
  dat_outhospital_diagnose <- na.omit(dat_outhospital_diagnose) %>% as.character()
  dat_dismiss_diagnose <- ifelse(str_detect(dat_outhospital_diagnose, "重度抑郁"),"重度抑郁",
                                 ifelse(str_detect(dat_outhospital_diagnose, "复发性抑郁"), "复发性抑郁",
                                        ifelse(str_detect(dat_outhospital_diagnose, "心境障碍"), "心境障碍",
                                               ifelse(str_detect(dat_outhospital_diagnose, "抑郁发作"), "抑郁发作",'其他'))))
  dat_dismiss_psycho <- ifelse(str_detect(dat_outhospital_diagnose, "不伴精神病"),'No',
                               ifelse(str_detect(dat_outhospital_diagnose, "，伴精神病"), 'Yes',
                                      NA))
  dat_dismiss_soma <- ifelse(str_detect(dat_outhospital_diagnose, "不伴躯体化"),'No',
                             ifelse(str_detect(dat_outhospital_diagnose, "，伴躯体化"), 'Yes',
                                    ifelse(str_detect(dat_outhospital_diagnose, "躯体症状障碍"), 'Yes',
                                           NA)))
  
  testPDF_2 <- str_replace_all(testPDF, pattern = '\n', replacement = '')
  testPDF_2 <- paste(testPDF_2, collapse = " ")
  
  ## enroll frequency
  dat_hosp_freq <- str_extract(testPDF_2, pattern = "第.*次住院")[1]
  dat_hosp_freq2 <- str_extract(dat_hosp_freq, pattern = "[0-9]{1,3}")
  
  ## in hospital duration
  dat_hpospital_day <- str_extract(testPDF_2, pattern = "住院天数：.*天科")[1]
  dat_hpospital_day2 <- str_extract(dat_hpospital_day, pattern = "[0-9]{1,3}")
  
  dat_suicide_risk <- str_extract(testPDF_2, pattern = "自杀风险评估为.风险")
  dat_suicide_risk <- na.omit(dat_suicide_risk) %>% as.character()
  if (is.empty(dat_suicide_risk)) {
    dat_suicide_risk <- str_extract(testPDF_2, pattern = "自杀风险.风险")
    dat_suicide_risk <- na.omit(dat_suicide_risk) %>% as.character()
  }
  if (is.empty(dat_suicide_risk)) {
    dat_suicide_risk <- str_extract(testPDF_2, pattern = "自杀风险为.风险")
    dat_suicide_risk <- na.omit(dat_suicide_risk) %>% as.character()
  }
  if (is.empty(dat_suicide_risk)) {
    dat_suicide_risk <- str_extract(testPDF_2, pattern = "评估为.风险")
    dat_suicide_risk <- na.omit(dat_suicide_risk) %>% as.character()
  }
  dat_suicide_risk_rank <- ifelse(str_detect(dat_suicide_risk, pattern = "低风险"), 'low',
                                  ifelse(str_detect(dat_suicide_risk, pattern = "中风险"), 'medium',
                                         'high'))
  
  ## collect the information
  
  dat_table <- data.frame(
    name = dat_name2,
    gender = dat_gender2,
    hospitalID = dat_hospID,
    hospitalFreq = dat_hosp_freq2,
    Date_Enroll = dat_enroll_date2,
    Date_Dismiss = dat_dismiss_date2,
    hospital_duration = dat_hpospital_day2,
    enroll_depression = dat_enroll_diagnose,
    enroll_psycho = dat_enroll_psycho,
    enroll_somatization = dat_enroll_soma,
    dismiss_depression = dat_dismiss_diagnose,
    dismiss_psycho = dat_dismiss_psycho,
    dismiss_somatization = dat_dismiss_soma,
    dismiss_suicide_risk = dat_suicide_risk_rank
  )
  
  return(dat_table)
}


ExtractEnrollInformation_reshape <- function(path_to_pdf){
  library(pdftools)
  library(tidyverse)
  
  # load the pdf 
  testPDF <- pdf_text(path_to_pdf)
  
  # subject name
  dat_name <- str_extract(testPDF[1], pattern = "姓名：.{2,3}")
  dat_name2 <- str_split(dat_name, pattern = '：') %>% unlist() %>% .[2] %>% 
    trimws(which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  # job/occupation
  dat_job <- str_extract(testPDF[1], pattern = "职业：.{2,8}")
  dat_job2 <-  dat_job %>% str_split(pattern = '：') %>% unlist() %>% .[2] %>% 
    trimws(which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  # marriage
  dat_marriage <- str_extract(testPDF[1], pattern = "婚姻状况：.{2}")
  dat_marriage2 <- dat_marriage %>% str_split(pattern = '：') %>% unlist() %>% .[2] %>% 
    trimws(which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  # nationality 
  dat_nation <- str_extract(testPDF[1], pattern = "民族：.{2}")
  dat_nation2 <- dat_nation %>% str_split(pattern = '：') %>% unlist() %>% .[2] %>% 
    trimws(which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  # education 
  dat_education <- str_extract(testPDF[1], pattern = "文化程度：.{2}")
  dat_education2 <- dat_education %>% str_split(pattern = '：') %>% unlist() %>% .[2] %>% 
    trimws(which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  # birthplace
  dat_birthplace <- str_extract(testPDF[1], pattern = "生地：.{7}") 
  dat_birthplace2 <- dat_birthplace %>% str_split(pattern = '：') %>% unlist() %>% .[2] %>% 
    trimws(which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  # self narrative
  testPDF2 <- str_replace_all(testPDF, pattern = '\n', replacement = '')
  dat_selfrep <- str_extract(testPDF2, pattern = "主诉：.*现病史：")
  dat_selfrep <- na.omit(dat_selfrep) %>% as.character()
  
  # current illness
  dat_curIllness <- str_extract(testPDF2, pattern = "现病史：.* ") 
  dat_curIllness <- na.omit(dat_curIllness) %>% as.character()
  
  # past history
  dat_pastHist <- str_extract(testPDF2, pattern = "既往史：.*个人史")
  dat_pastHist <- na.omit(dat_pastHist) %>% as.character()
  
  # individual history
  dat_indiviHist <- str_extract(testPDF2, pattern = "个人史.* ")
  dat_indiviHist <- na.omit(dat_indiviHist) %>% as.character()
  
  dat_table <- data.frame(
    name = dat_name2,
    occupation = dat_job2,
    education = dat_education2,
    nationality = dat_nation2,
    marriage = dat_marriage2,
    birthplace = dat_birthplace2
  )
  
  return(dat_table)
}


ExtractDismissInformation_raw <- function(path_to_pdf){
  library(pdftools)
  library(tidyverse)
  
  # load the pdf 
  testPDF <- pdf_text(path_to_pdf)
  
  ## subj name
  dat_name <- str_extract(testPDF, pattern = "姓名：.{2,3}")[1]
  dat_name2 <- str_split(dat_name, pattern = '：') %>% unlist() %>% .[2] %>% 
    trimws(which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
    
  ## hospital ID
  dat_hospID <- str_extract(testPDF, pattern = "[0-9]{10}")[1]
  
  ## gender
  dat_gender <- str_extract(testPDF, pattern = "[男,女]")[1]
  
  ## enroll time
  dat_enroll_date <- str_extract(testPDF, pattern = "入院时间：[0-9]{4}-[0-9]{2}-[0-9]{2}")[1]
  
  ## dismiss time
  dat_dismiss_date <- str_extract(testPDF, pattern = "出院时间：[0-9]{4}-[0-9]{2}-[0-9]{2}")[1]

    ## enroll diagnose
  dat_inhospital_diagnose <- str_extract(testPDF, pattern = "入院诊断：.*\n")[1]
  
  testPDF <- paste0(testPDF)

  ## dismiss diagnose
  dat_outhospital_diagnose <- str_extract(testPDF, pattern = "出院诊断：.*\n")
  dat_outhospital_diagnose <- na.omit(dat_outhospital_diagnose) %>% as.character()
  
  testPDF_2 <- str_replace_all(testPDF, pattern = '\n', replacement = '')
  testPDF_2 <- paste(testPDF_2, collapse = " ")
  
  ## enroll frequency
  dat_hosp_freq <- str_extract(testPDF_2, pattern = "第.*次住院")[1]
  
  ## in hospital duration
  dat_hpospital_day <- str_extract(testPDF_2, pattern = "住院天数：.*天科")[1]
  
  ## dismiss suicide risk aessessment
  dat_suicide_risk <- str_extract(testPDF_2, pattern = "自杀风险评估为.风险")
  dat_suicide_risk <- na.omit(dat_suicide_risk) %>% as.character()
  if (is.empty(dat_suicide_risk)) {
    dat_suicide_risk <- str_extract(testPDF_2, pattern = "自杀风险.风险")
    dat_suicide_risk <- na.omit(dat_suicide_risk) %>% as.character()
  }
  if (is.empty(dat_suicide_risk)) {
    dat_suicide_risk <- str_extract(testPDF_2, pattern = "自杀风险为.风险")
    dat_suicide_risk <- na.omit(dat_suicide_risk) %>% as.character()
  }
  if (is.empty(dat_suicide_risk)) {
    dat_suicide_risk <- str_extract(testPDF_2, pattern = "评估为.风险")
    dat_suicide_risk <- na.omit(dat_suicide_risk) %>% as.character()
  }

  ## collect the information
  
  dat_table <- data.frame(
    name = dat_name2,
    gender = dat_gender,
    hospitalID = dat_hospID,
    hospitalFreq = dat_hosp_freq,
    Date_Enroll = dat_enroll_date,
    Date_Dismiss = dat_dismiss_date,
    hospital_duration = dat_hpospital_day,
    enroll_diagnose = dat_inhospital_diagnose,
    dismiss_diagnose = dat_outhospital_diagnose,
    dismiss_suicide_risk = dat_suicide_risk
  )
  
  return(dat_table)
}


ExtractEnrollInformation_raw <- function(path_to_pdf){
  library(pdftools)
  library(tidyverse)
  
  # load the pdf 
  testPDF <- pdf_text(path_to_pdf)
  
  # subject name
  dat_name <- str_extract(testPDF[1], pattern = "姓名：.{2,3}")
  
  # job/occupation
  dat_job <- str_extract(testPDF[1], pattern = "职业：.{2,8}")
  
  # marriage
  dat_marriage <- str_extract(testPDF[1], pattern = "婚姻状况：.{2}")
  
  # nationality 
  dat_nation <- str_extract(testPDF[1], pattern = "民族：.{2}")
  
  # education 
  dat_education <- str_extract(testPDF[1], pattern = "文化程度：.{2}")
  
  # birthplace
  dat_birthplace <- str_extract(testPDF[1], pattern = "生地：.{7}") 
  
  # self narrative
  testPDF2 <- str_replace_all(testPDF, pattern = '\n', replacement = '')
  testPDF2 <- paste(testPDF2, collapse = " ")
  
  dat_selfrep <- str_extract(testPDF2, pattern = "主诉：.*现病史：")
  dat_selfrep <- na.omit(dat_selfrep) %>% as.character()
  
  # current illness
  dat_curIllness <- str_extract(testPDF2, pattern = "现病史：.*既往史：") 
  dat_curIllness <- na.omit(dat_curIllness) %>% as.character()
  
  # past history
  dat_pastHist <- str_extract(testPDF2, pattern = "既往史：.*个人史")
  dat_pastHist <- na.omit(dat_pastHist) %>% as.character()
  
  # individual history
  dat_indiviHist <- str_extract(testPDF2, pattern = "个人史.* ")
  dat_indiviHist <- na.omit(dat_indiviHist) %>% as.character()
  
  dat_table <- data.frame(
    name = dat_name,
    occupation = dat_job,
    education = dat_education,
    nationality = dat_nation,
    marriage = dat_marriage,
    birthplace = dat_birthplace,
    self_report = dat_selfrep,
    current_illness = dat_curIllness,
    past_illness = dat_pastHist,
    autobiography = dat_indiviHist
  )
  
  return(dat_table)
}

