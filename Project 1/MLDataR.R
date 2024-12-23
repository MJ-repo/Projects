
rm(list=ls())
install.packages("MLDataR")
install.packages("caTools")
library(MLDataR)
library(dplyr)
library(tidyr)
library(tidymodels)
library(data.table)
library(ConfusionTableR)
library(OddsPlotty)
library(caTools)
glimpse(MLDataR::heartdisease)
head(heartdisease)
str(heartdisease)

dat<-as.data.frame(heartdisease)
head(dat[2],10)
dat$Sex<-as.factor(dat$Sex)
str(dat)

hd <- heartdisease %>%
  mutate(across(where(is.character), as.factor),
         HeartDisease = as.factor(HeartDisease)) |>
  na.omit()
str(hd)


set.seed(123)
split_prop <- 0.8
testing_prop <- 1 - split_prop
split <- rsample::initial_split(hd, prop = split_prop)
training <- rsample::training(split)
testing <- rsample::testing(split)
# Print a custom message to show the samples involved
training_message <- function() {
  message(
    cat(
      'The training set has: ',
      nrow(training),
      ' examples and the testing set has:',
      nrow(testing),
      '.\nThis split has ',
      paste0(format(100 * split_prop), '%'),
      ' in the training set and ',
      paste0(format(100 * testing_prop), '%'),
      ' in the testing set.',
      sep = ''
    )
  )
}
training_message()
## The training set has: 734 examples and the testing set has:184.
## This split has 80% in the training set and 20% in the testing set.
## 

lr_hd_fit <- logistic_reg() %>%
  set_engine("glm") %>% 
  set_mode("classification") %>% 
  fit(HeartDisease ~ ., data = training)

tidy(lr_hd_fit)
tidy(lr_hd_fit) %>% 
  filter(p.value < 0.05) %>% 
  pull(term)

tidy_oddsplot <- OddsPlotty::odds_plot(
  lr_hd_fit$fit,
  title = "Heart Disease Odds Plot",
  point_col = "#6b95ff",
  h_line_color = "red"
)
tidy_oddsplot <- tidy_oddsplot$odds_plot +
  theme(legend.position = "none") +
  geom_text(
    label = round(tidy_oddsplot$odds_plot$data$OR, digits = 3),
    hjust = -0.5,
    vjust = 1,
    cex = 2.8
  )
tidy_oddsplot

library(ConfusionTableR)
# Use our model to predict labels on to testing set
predictions <- cbind(predict(lr_hd_fit, new_data = testing),
                     testing)
# Create confusion matrix and output to record level for storage to monitor concept drift
cm <- ConfusionTableR::binary_class_cm(
  predictions$.pred_class,
  predictions$HeartDisease,
  mode = 'everything',
  positive = '1'
)
# Access the confusion matrix list object
cm$confusion_matrix
