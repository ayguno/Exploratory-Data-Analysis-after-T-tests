# Exploratory-Data-Analysis-after-T-tests
Facilitates the post processing of the data table obtained from the Shiny server ModT test for a TMT10 experiment

It takes the shiny server result table and class vector BOTH CONVERTED TO CSV files as input (char) 

Project = give a specific name for the project (char)

Pvalue = the adj.p value cut off specified to determine significance (num)

xlimit, ylimit = desired magnitude boundaries in the scatterplots (num)

prints various interesting tables including what proteins deemed to be significantly up/down regulated (adj.p value <Pv) in all of the 4 test conditions

It also marks the significant hits on the log2 scatterplots and prints the gene symbol of the top 10 (based on the adj p-value) data points

It also prepares two files to facilitate the GSEA (java and R-based (ssGSEA), respectively)

Last Update on 09/13/2016 : added, check.names=F in the gct preparation to avoid additional letters,
also added Rep1 and Rep2 strings at the end of the column names as they are generated to ensure the variable uniqueness
