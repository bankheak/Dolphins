library(vegan)
library(tidyverse)

help(vegdist)
data(varespec)
str(varespec)
help("varespec")
t(varespec)



# select variables from the raw data
aux = orig_data[1:1000, 
                c('Code', 'Behaviors', 'HumanInteraction', 'ConfHI')]

# Use 'Behaviors' variable to extract "Feed" and create another variable with two classes (Feed, Other)

aux$Foraging = "Other"
aux$Foraging[grepl(pattern = 'Feed',
                    x = aux$Behaviors,
                    ignore.case = FALSE, perl = FALSE,
                    fixed = FALSE, useBytes = FALSE)] = "Feed"


IDbehav = table(aux$Code, aux$Foraging)

rawHI = as.matrix(table(aux$Code, aux$ConfHI)[,2:6])

IDdata = data.frame(Foraging = IDbehav[,1])
IDdata$HI = as.vector(rowSums(rawHI))
IDdata$HIprop = IDdata$HI/IDdata$Foraging


IDdata = as.data.frame(table(aux$Code))
IDdata$HI = as.vector(rowSums(rawHI))
IDdata$HIprop = IDdata$HI/IDdata$Freq

identical(row.names(rawHI), IDdata$Var1)


# Count how many times each individual was 'Feed'
# Count how many times each ind was in HI or not
# Divide the two to create the 'propHI'


aux$Behaviors[grepl(pattern = 'Feed',
                    x = aux$Behaviors,
                    ignore.case = FALSE, perl = FALSE,
                    fixed = FALSE, useBytes = FALSE)]

as.vector(table(aux$Behaviors))
table(aux$HumanInteraction)


# pick the feeding events
aux[which(aux$Behaviors == 'pFeed'), ]

# 
aux[which(aux$HumanInteraction != 999), ]




# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
fake_HIprop = t(varespec)[,1]
as.matrix(vegdist(fake_HIprop, method = 'euclidean'))
