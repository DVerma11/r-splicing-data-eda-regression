################### PROJECT ASSSIGNMENT####################
library(ggplot2)
############################## 1. DATA EXPLORATION ##############################
#1.1 Load the provided dataset, "splicing_data.csv," into R- no need to convert in long data- linear reg. based on slicing factor
splicingData <- read.csv(file.choose(), header= TRUE)

#1.2 Perform summary statistics of the dataset to get an overview.
str(splicingData)
summary(splicingData)

#1.3 Missing values and outliers
#A. Check for missing values 
sum(is.na(splicingData))

#B Check for outliers: using boxplot
splicingData_Boxplot <-boxplot(splicingData[2:5])
print(splicingData_Boxplot)

#1.4 Generate a correlation matrix to assess the relationships between the splicing factors and the splicing event.

#1.4.1 extract numeric variables
splicingDataMatrix <- as.matrix(splicingData[, c("SplicingFactor1", "SplicingFactor2", "SplicingFactor3","SplicingEvent" )])

#1.4.2 Using rcorr()- default is pearson
library(Hmisc)
rcorr(splicingDataMatrix)

########################## 2. DATA VISUALIZATION ##############################

#2.1 Create visualizations to explore the distribution of each splicing factor's expression and the splicing event.

## Using Histogram for each variable with normal curve
#SplicingFactor1  
SplicingFactor1.hist <- ggplot(splicingData, aes(x=SplicingFactor1)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  labs(title="Graph 1: Distribution of Splicing Factor 1", x="SPLICING FACTOR 1", y="Density",
       caption = "Image by Divya Verma")
SplicingFactor1.hist  + 
  stat_function(fun = dnorm, 
                args = list(mean = mean(splicingData$SplicingFactor1, na.rm = TRUE), 
                            sd = sd(splicingData$SplicingFactor1, na.rm = TRUE)), 
                colour = "black", 
                size = 1)

#SplicingFactor2 
SplicingFactor2.hist <- ggplot(splicingData, aes(x=SplicingFactor2)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  labs(title="Graph 2: Distribution of Splicing Factor 1", x="Splicing Factor 2", y="Density",
       caption = "Image by Divya Verma")
SplicingFactor2.hist  + 
  stat_function(fun = dnorm, 
                args = list(mean = mean(splicingData$SplicingFactor2, na.rm = TRUE), 
                            sd = sd(splicingData$SplicingFactor2, na.rm = TRUE)), 
                colour = "black", 
                size = 1)

#SplicingFactor3 
SplicingFactor3.hist <- ggplot(splicingData, aes(x=SplicingFactor3)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  labs(title="Graph 3: Distribution of Splicing Factor 3", x="Splicing Factor 3", y="Density",
       caption = "Image by Divya Verma")
SplicingFactor3.hist  + 
  stat_function(fun = dnorm, 
                args = list(mean = mean(splicingData$SplicingFactor3, na.rm = TRUE), 
                            sd = sd(splicingData$SplicingFactor3, na.rm = TRUE)), 
                colour = "black", 
                size = 1)

#SplicingEvent 
SplicingEvent.hist <- ggplot(splicingData, aes(x=SplicingEvent)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  labs(title="Graph 4: Distribution of Splicing Event", x="Splicing Event", y="Density",
       caption = "Image by Divya Verma")
SplicingEvent.hist  + 
  stat_function(fun = dnorm, 
                args = list(mean = mean(splicingData$SplicingEvent, na.rm = TRUE), 
                            sd = sd(splicingData$SplicingEvent, na.rm = TRUE)), 
                colour = "black", 
                size = 1)

#Quantifying Normality with Numbers(Assessing Skew and Kurtosis)
library(pastecs )
stat.desc(cbind(splicingData$SplicingFactor1, splicingData$SplicingFactor2, 
                splicingData$SplicingFactor3, splicingData$SplicingEvent), 
          basic = FALSE, norm = TRUE
) 
#Assessing Normality using Shapiro-Wilk normality test
shapiro.test(splicingData$SplicingFactor1) # p-value = 0.9349 [p>0.05]: non significantly/normal
shapiro.test(splicingData$SplicingFactor2) # p-value = 0.03691 [p<0.05]: significant/ non-normal
shapiro.test(splicingData$SplicingFactor3) # p-value = 0.06513 [p>0.05]: non significantly/normal
shapiro.test(splicingData$SplicingEvent) # p-value = 0.3462 [p>0.05]: non significant/normal


#2.2 Create graphs to visualize the relationships between each individual splicing factor and the splicing event.
## scatterplots to visualize relationship between variables- Splicing factor and splicing event.
#SplicingFactor1
scatter_SplicingFactor1<- ggplot(splicingData, aes(SplicingFactor1, SplicingEvent))
scatter_SplicingFactor1 + 
  geom_point() + 
  geom_smooth(method= "lm", colour= "blue") +
  labs(x="Splicing Factor 1", y= "Splicing Event", 
       title = "Graph 5: Relationship between Splicing Factor 1 and Splicing Event",
       caption="Image by Divya Verma")
#SplicingFactor2
scatter_SplicingFactor2<- ggplot(splicingData, aes(SplicingFactor2, SplicingEvent))
scatter_SplicingFactor2 + 
  geom_point() + 
  geom_smooth(method= "lm", colour= "red") +
  labs(x="Splicing Factor 2", y= "Splicing Event", 
       title = "Graph 6: Relationship between Splicing Factor 2 and Splicing Event",
       caption="Image by Divya Verma")
#SplicingFactor3
scatter_SplicingFactor3<- ggplot(splicingData, aes(SplicingFactor3, SplicingEvent))
scatter_SplicingFactor3 + 
  geom_point() + 
  geom_smooth(method= "lm", colour= "green") +
  labs(x="Splicing Factor 3", y= "Splicing Event", 
       title = "Graph 7: Relationship between Splicing Factor 3 and Splicing Event",
       caption="Image by Divya Verma")
#########################################
#3. PREDICTIVE MODELING
#3.1 Split the dataset into a training set (70%) and a testing set (30%).
#3.1.1 Create object n which is no. of rows of splicing data(to divide the data)
n <- nrow(splicingData)
n
#3.1.2 Divide into 70 %
train_indices <- sample(1:nrow(splicingData), 0.7 * nrow(splicingData))

#3.1.3 Create train data 
train_data <- splicingData[train_indices, ]

#3.1.4 Create test data 
test_data <- splicingData[-train_indices, ]

# 3.2 Build predictive Training model (linear regression)
TrainingModel <- lm(SplicingEvent ~ SplicingFactor1 + SplicingFactor2 + SplicingFactor3, 
                      data=train_data)
summary(TrainingModel) 

#3.3 Evaluate the model's performance using appropriate metrics on the testing set.
#3.3.1 Evaluate model performance on the testing set
predictedTest <- predict(TrainingModel, newdata = test_data)
#3.3.2 Calculate R-squared  
R_squaredTest <- cor(predictedTest, test_data$SplicingEvent)^2
R_squaredTest
#Result: 0.6954288

# 3.3.3 Calculate Mean Squared Error (MSE) 
MSE_Test <- mean((predictedTest - test_data$SplicingEvent)^2)
MSE_Test
# Result: 3.105155

#3.4 Interpret the model's coefficients 

coefficients_TrainingModel <- coef(TrainingModel)
coefficients_TrainingModel

