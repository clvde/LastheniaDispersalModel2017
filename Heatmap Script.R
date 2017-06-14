library(ggplot2)

#------------------
# CREATE DATA FRAME
#------------------
df.team_data <- expand.grid(teams = c("Team A", "Team B", "Team C", "Team D"), metrics = c("Metric 1", "Metric 2", "Metric 3", "Metric 4", "Metric 5"))
df.team_data
# add variable: performance
set.seed(41)
df.team_data$performance <- rnorm(nrow(df.team_data))

#inspect
head(df.team_data)

#---------------------------
# PLOT: heatmap
# - here, we use geom_tile()
#---------------------------

ggplot(data = df.team_data, aes(x = metrics, y = teams)) + geom_tile(aes(fill = performance)) 

#---------------------------
# CREATE DATA FRAME
#---------------------------


