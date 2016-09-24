
#This function is aimed to aid the post processing of the data table obtained from the Shiny server ModT test for a TMT10 experiment
# It takes the shiny server result table and class vector BOTH CONVERTED TO CSV files as input (char) 
# Project = give a specific name for the project (char)
# Pvalue = the adj.p value cut off specified to determine significance (num)
# xlimit, ylimit = desired magnitude boundaries in the scatterplots (num)
# and prints various interesting tables including what proteins deemed to be
# significantly up/down regulated (adj.p value <Pv) in all of the 4 test conditions
# It also marks the significant hits on the log2 scatterplots and prints the gene symbol of the top 10 (based on the adj p-value) data points
# It also prepares two files to facilitate the GSEA (java and R-based (ssGSEA), respectively)

#Last Update on 09/13/2016 : added, check.names=F in the gct preparation to avoid additional letters,
# also added Rep1 and Rep2 strings at the end of the column names as they are generated to ensure the variable uniqueness

postModT <-function(ModT,ClassV,Pvalue,xlimit,ylimit,Project) {
  
  ## reading ModT, Class Vector tables, specified adj p value cut off and Project identifier
  
  MT<-read.csv(ModT,header=TRUE)
  ClV<-read.csv(ClassV,header=TRUE)
  Pv<-Pvalue
  Pr<-Project
  xl<-xlimit
  yl<-ylimit
  
  ## labeling the class vector with experiment names
  
  ClV<-data.frame(ClV, Repli=c("Rep1","Rep2","Rep1","Rep2","Rep1","Rep2","Rep1","Rep2"))
  
  
  ##Printing the significantly regulated proteins (adj.p. val<Pv) for each experiment
  
  print("Printing the significantly regulated proteins (adj.p. val<Pv) for each experiment")
  
  #Looping over the individual experiments 
  
  Eindx<-NULL
  SumE <-NULL
  SumEUP<-NULL
  SumEDN<-NULL
  
  a<-1
  n<-nrow(ClV)/2
  
  for(i in 1: n) {
    
    
    adjp<-MT[,paste("adj.P.Val", ClV[a,2],sep=".")]
    logFC<-MT[,paste("logFC", ClV[a,2],sep=".")]
    
    E<-which(adjp<Pv)
    EUP<-which(adjp<Pv & logFC>0) 
    EDN<-which(adjp<Pv & logFC<0)
    
    Exp1<- MT[E,]
      write.csv(file=paste(Pr,ClV[a,2],"adjPVal",Pv, "significant",".csv",sep="_"),Exp1)
    Exp1UP<-MT[EUP,]
      write.csv(file=paste(Pr,ClV[a,2],"adjPVal",Pv,"significant","UP",".csv",sep="_"),Exp1UP)
    Exp1DN<-MT[EDN,]
      write.csv(file=paste(Pr,ClV[a,2],"adjPVal",Pv,"significant","DOWN",".csv",sep="_"),Exp1DN)
    
    #Summarizing the number of significantly regulated proteins for each experiment as the loop moves
      
      Eindx<-c(Eindx,as.character(ClV[a,2]))
      SumE <-c(SumE,length(E))
      SumEUP<-c(SumEUP,length(EUP))
      SumEDN<-c(SumEDN,length(EDN))
     
    a<-a+2
  }
  
  
  signsum<-data.frame(Eindx,SumE,SumEUP,SumEDN)
  colnames(signsum)<-c("Experiment",paste("Total # of significantly regulated proteins","adj.p.val<",Pv,Sep=""), "UP", "DOWN")
  write.csv(file=paste(Pr,"summary","adjPVal",Pv,"numberofsign_reg_proteins.csv",sep="_"),signsum)
 
  print("The significantly regulated proteins sucessfully printed")
  
  
  ##Next, preparing scatterplots where the significant hits are color marked
  
  n<-nrow(ClV)/4
  
  par(mfrow=c(2,n), mar=c(5,5,5,5)) # plotting  scatter plots for the  TMT10 ratio pairs
  
  a<-1
  n1<-nrow(ClV)/2
  
  for(i in 1:n1){
    
    adjp<-MT[,paste("adj.P.Val", ClV[a,2],sep=".")]
    
    x<-which(colnames(MT)==ClV[a,1])
    y<-which(colnames(MT)==ClV[a+1,1])
    
    plot(MT[,x],MT[,y],pch=20, col=ifelse(adjp<Pv,"red","black"), cex=1,xlim=c(-xl,xl),ylim=c(-yl,yl),main=ClV[a,2],cex.lab=0.4, xlab=colnames(MT)[x],ylab=colnames(MT)[y])
    abline(v=0,h=0,lwd=1,lty=2)
                                  
    a<-a+2
  }
  
        pdf(paste(Pr,"adjPVal",Pv,"Significants_are_red_scatterplots.pdf",sep="_")) # Printing the PDF version of the scatterplots
        
        n<-nrow(ClV)/4
        par(mfrow=c(2,n), mar=c(5,5,5,5)) # plotting  scatter plots for the  TMT10 ratio pairs
        
        a<-1
        n1<-nrow(ClV)/2
        
        for(i in 1:n1){
          
          adjp<-MT[,paste("adj.P.Val", ClV[a,2],sep=".")]
          
          x<-which(colnames(MT)==ClV[a,1])
          y<-which(colnames(MT)==ClV[a+1,1])
          
          plot(MT[,x],MT[,y],pch=20, col=ifelse(adjp<Pv,"red","black"), cex=1,xlim=c(-xl,xl),ylim=c(-yl,yl),main=ClV[a,2],cex.lab=0.4, xlab=colnames(MT)[x],ylab=colnames(MT)[y])
          abline(v=0,h=0,lwd=1,lty=2)
          text(paste("adj.p.val<",Pv,sep=""),x=1,y=-1.8,col="Red")
          
          a<-a+2
        }
        
        dev.off()
        
        print("The significantly regulated proteins sucessfully plotted")
        
        print("Labeling the names of the top 10 significantly regulated proteins on the plots")
        
  ##Next, marking the gene symbols of the top 10 up/down regulated proteins on scatterplots where the significant hits are color marked
    
     
           
    pdf(file=paste(Pr,"Top10geneslabelledsignificantredscatterplots.pdf",sep="_")) # these 4 plots will be printed into PDF
    
   
    require("maptools") #To use pointlabel function below, very useful in avoiding overlapping data labels 
    
    a<-1
    n1<-nrow(ClV)/2
    
    for(i in 1:n1){
      
      adjp<-MT[,paste("adj.P.Val", ClV[a,2],sep=".")]
      logFC<-MT[,paste("logFC", ClV[a,2],sep=".")]
      
      x<-which(colnames(MT)==ClV[a,1])
      y<-which(colnames(MT)==ClV[a+1,1])
      
      plot(MT[,x],MT[,y],pch=20, col=ifelse(adjp<Pv,"red","black"), cex=1,xlim=c(-xl,xl),ylim=c(-yl,yl),main=ClV[a,2],cex.lab=0.4, xlab=colnames(MT)[x],ylab=colnames(MT)[y])
      abline(v=0,h=0,lwd=1,lty=2) 
      text(paste("adj.p.val<",Pv,sep=""),x=1,y=-1.8,col="Red")
     
      
            w<-which(adjp<Pv) #For the first experiment in the loop, find the significant rows
            temp<-MT[w,] #Collect the significant values for this experiment into a new temporary data frame
            tlogFC<-which(colnames(temp)== paste("logFC", ClV[a,2],sep="."))
            sorttemp<-temp[order(temp[,tlogFC]),] # Ascending sort the 'temp' based on LogFC of the specific experiment
            
            tUP<-which(sorttemp[,tlogFC]>0) #Extract the upregulated fraction of the sorted 'temp'
            sorttempUP<-sorttemp[tUP,]
            
            tDN<-which(sorttemp[,tlogFC]<0)#Extract the downregulated fraction of the sorted 'temp'
            sorttempDN<-sorttemp[tDN,]
            
            tempUP<-tail(sorttempUP,10) #Extract the 10 most significantly upregulated proteins
            if(nrow(tempUP)>0) 
              pointLabel(tempUP[,x],tempUP[,y],labels=tempUP[,"geneSymbol"],cex=0.75) #Label them in the existing plot
          
            tempDN<-head(sorttempDN,10) #Extract the 10 most significantly downregulated proteins
            if(nrow(tempDN)>0)
              pointLabel(tempDN[,x],tempDN[,y],labels=tempDN[,"geneSymbol"],cex=0.75) #Label them in the existing plot
            
      a<-a+2
    } 
    
    dev.off()
    
##Next, generating the files required for the GSEA analyses
    
    print("Generating the files required for the GSEA analyses")
    
    #Generating the .gct file required for the ssGSEA analysis in R environment
    
    gct<-data.frame(NAME=MT[,"geneSymbol"],Description=MT[,"id"],check.names = F) # First collecting Uniprot and Gene symbol columns from MT table
    
    a<-1
    n1<-nrow(ClV)/2
    
    P1<-which(colnames(MT)== paste("P.Value", ClV[1,2],sep="."))
    PL<-which(colnames(MT)== paste("P.Value", ClV[nrow(ClV),2],sep="."))
    
    for(i in 1:n1){
      
      x<-which(colnames(MT)==ClV[a,1])
      y<-which(colnames(MT)==ClV[a+1,1])
      
          gct<-data.frame(gct,MT[,x],MT[,y],check.names = F) #Appending the TMT Log2 ratios for two replicates into gct
          colnames(gct)[a+2]<-paste(ClV[a,2],ClV[a,3],"Rep1",sep="_") #Column name of Rep1
          colnames(gct)[a+3]<-paste(ClV[a+1,2],ClV[a+1,3],"Rep2",sep="_") #Column name of Rep2
          
          a<-a+2 }
    
    ### Unique gene name filter: next filter any repeated gene names from the list based on nominal p-values: keep the gene name that has the lowest average p-value among the repeated gene names  
    
    AvP<-apply(MT[,P1:PL],1,mean) # Get the average nominal p-value for each protein in the table
    
    gct[,"Av. P value"]<-AvP # Add the average nominal p-value column for each protein to the gct 
    
    sort_gct<-gct[order(gct[,"Av. P value"]),] # Ascending sort the 'gct' based on 'AvP' 
    
    ##Next filter the sorted gct based on gene names
    
    unique_gct<-sort_gct[!duplicated(sort_gct$NAME),] #Choosing the unique gene name which has the smallest average p-value among the duplicates
    
    as.character(unique_gct[,1])
    w<- which(unique_gct[,1]=="") # Find any genes with missing names
    
    unique_gct<-unique_gct[-w,] # Remove the rows of any genes with missing names
    
    unique_gct<-unique_gct[,-(ncol(unique_gct))] # Remove the Av. P value column we don't need in the final table

    l<-nrow(unique_gct) # get the number of rows in the unique_gct
    
    f<-file(paste(Pr,".gct",".csv",sep="_"),"w") # Open a connection and write the header to it first, then write the data frame. Note that the custom project name is printed to the .gct file name

    writeLines("#1.2",f) # Write the first generic line to the file
    
    writeLines(paste(l,8,sep=",\ "),f) # Write the second generic line to the file, which consists of # of rows in the data set and number of conditions that are individually tested in ssGSEA (by default, this number is 8 for a TMT10 experiment)
    
    write.csv(unique_gct,f) # Writes the 'unique_gct' with the lines above added to f, by using the connection f
    
    close(f) #Close the connection f. The desired ssGSEA file is ready, can be used after deleting the .csv in the data file
    
   
    
     ##Next prepare the GSEA files required for the java-based GSEA analysis
    
    
    #First generating the ratio table in delimited .txt format
    
    lineargct<-2^(unique_gct[,3:ncol(unique_gct)]) #linearized the log2  ratios
    
    lineargct<-data.frame(geneSymbol=unique_gct[,1],lineargct) #add the gene name column
    
    #Finally add two additional columns as relative change control, all containing value of 1
    
    control_rep<-seq(1,1,length.out = l)
    
    lineargct<-data.frame(lineargct,control=control_rep,control=control_rep)
    
        #Writing the ratio table for java-based ssGSEA analysis. Convert this to txt file before use
        write.csv(file=paste(Pr,"linearized_ratio_table_for_Java_GSEA",".txt",".csv",sep="_"),lineargct)
        
    #Prepare and print the .cls table required for the ava-based ssGSEA analysis  
    
    df<-data.frame(NULL)
        
    f<-file(paste(Pr,"cls_table_for_Java_GSEA",".cls",".csv",sep="_"),"w") # Open a connection
    
    writeLines(paste(ncol(lineargct)-1,(ncol(lineargct)-1)/2,1,sep=",\ "),f) # Write the first generic line to the file
    
    temp<-paste("#",ClV[1,2],sep="")
    
    writeLines(paste(temp,ClV[3,2],ClV[5,2],ClV[7,2],"Control",sep=",\ "), f) #Write the second generic line to the file

    writeLines(paste(ClV[1,2],ClV[2,2],ClV[3,2],ClV[4,2],ClV[5,2],ClV[6,2],ClV[7,2],ClV[8,2],"Control","Control",sep=",\ "), f) #Write the third generic line to the file

    write.csv(df,f) # Writes the file
    
    close(f) #Close the connection f.
    
    dev.off()
    
    }  