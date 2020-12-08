###########################################################################################
########### Calculate mean roughness pointcloud from ClouCompare   #######################
########################## Update:  Dec.1, 2020 ###########################################
library(Morpho)
library(Arothron)
library(Rvcg)
library(rgl)
library(reshape2)
library(tidyverse)
library(stringr) # remove characters
library(dplyr)

#1. Extract only the outside of the oyster clusters.
#2. Calculate: surface area, roughness, 
#3. Calculate: difference between CT vs. Arothron1.0 vs. Arothron3.0

### 1. Extract only the outside of the oyster clusters

  setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh") # put your data file in this folder.or change directory

    (Full.cluster_mesh_folder <- list.files(pattern = ".ply"))
  Full.cluster_mesh_folder.select<-Full.cluster_mesh_folder[c(30)] # this is used to subset a part of data for pilot analyses

  Arothron.clusters.distance.POVnumber<-function(putanythinghere){
  start_time<-Sys.time()
  
        for (j in unique(Full.cluster_mesh_folder.select)){ 
          oyster.cluster<-vcgPlyRead(j)
          #
          # full.cluster<-ply2mesh(j) # only for # 4 cluster
          j.name<-unlist(strsplit(j,"_mesh.ply"))

        
        #  1.1 Arothron 
            calse_ectocast<-ext.int.mesh(mesh=oyster.cluster, views=160, dist.sphere=1, param1=3.0, default=TRUE,
                                           import_pov=NULL,matrix_pov=FALSE,expand=1, scale.factor=1, 
                                         num.cores=2,method="calse")
            vis_inv_ecto<-out.inn.mesh(calse_ectocast,oyster.cluster,plot=FALSE)
            vis_mesh<-vcgIsolated(vis_inv_ecto$visible)
            vis_mesh.clean<-vcgClean(vis_mesh,sel=c(0,1,2,3,4))
            
        # 1.2 Export file 
            setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist2")
            vcgPlyWrite(vis_mesh.clean,filename=paste0("A160_1_",j.name)) # export in ply or vcgObjWrite()
              print(paste0("A160_1_",j.name,"completed"))
              
              # 2.1 Arothron 
              calse_ectocast<-ext.int.mesh(mesh=oyster.cluster, views=160, dist.sphere=3, param1=3.0, default=TRUE,
                                           import_pov=NULL,matrix_pov=FALSE,expand=1, scale.factor=1, 
                                           num.cores=2,method="calse")
              vis_inv_ecto<-out.inn.mesh(calse_ectocast,oyster.cluster,plot=FALSE)
              vis_mesh<-vcgIsolated(vis_inv_ecto$visible)
              vis_mesh.clean<-vcgClean(vis_mesh,sel=c(0,1,2,3,4))
              
              # 2.2 Export file 
              setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist2")
              vcgPlyWrite(vis_mesh.clean,filename=paste0("A160_3_",j.name)) # export in ply or vcgObjWrite()
              print(paste0("A160_3_",j.name,"completed"))
        
        # 3.1 Arothron 
        calse_ectocast<-ext.int.mesh(mesh=oyster.cluster, views=80, dist.sphere=1, param1=3.0, default=TRUE,
                                     import_pov=NULL,matrix_pov=FALSE,expand=1, scale.factor=1, 
                                     num.cores=2,method="calse")
        vis_inv_ecto<-out.inn.mesh(calse_ectocast,oyster.cluster,plot=FALSE)
        vis_mesh<-vcgIsolated(vis_inv_ecto$visible)
        vis_mesh.clean<-vcgClean(vis_mesh,sel=c(0,1,2,3,4))
        
        # 3.2 Export file 
        setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist2")
        vcgPlyWrite(vis_mesh.clean,filename=paste0("A80_1_",j.name)) # export in ply or vcgObjWrite()
        print(paste0("A80_1_",j.name,"completed"))
              
              # 4.1 Arothron 
              calse_ectocast<-ext.int.mesh(mesh=oyster.cluster, views=80, dist.sphere=3, param1=3.0, default=TRUE,
                                           import_pov=NULL,matrix_pov=FALSE,expand=1, scale.factor=1, 
                                           num.cores=2,method="calse")
              vis_inv_ecto<-out.inn.mesh(calse_ectocast,oyster.cluster,plot=FALSE)
              vis_mesh<-vcgIsolated(vis_inv_ecto$visible)
              vis_mesh.clean<-vcgClean(vis_mesh,sel=c(0,1,2,3,4))
              
              # 4.2 Export file 
              setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist2")
              vcgPlyWrite(vis_mesh.clean,filename=paste0("A80_3_",j.name)) # export in ply or vcgObjWrite()
              setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh") # get back to the original folder destination
              print(paste0("A80_3_",j.name,"completed"))
    
        }
            end_time<-Sys.time()
            print(end_time-start_time) # estimate the time it took for the analyses
            
  }
  
  Arothron.clusters.distance.POVnumber(runit)
  
  
  
  

  
  
####################################################################################
#### 2. Calculate: distance from barycenter to the edge, surface area, and roughness


getwd()
setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist1")


(visib_folder <- list.files(pattern = ".ply"))
visib_folder.select<-visib_folder# [c(41:45,81:86,146:150,211:216)]

# make result saving matrix
R.result<-matrix(0,length(visib_folder.select),7)
rownames(R.result)<-strsplit(visib_folder.select,".ply");R.result
colnames(R.result)<-c("Cluster","POV","dist.sphere","Centroid.distance","SA","R.mean","R.SD")

# k<-"A160_3_Cluster39.ply"


    for (k in visib_folder.select){ # "A160_3_Cluster39.ply"){ # 
          #read ply file 
          setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist1")
          oyster.cluster<-vcgPlyRead(k)
         
          # full.cluster<-ply2mesh("Cluster4_mesh.ply") # only for # 4 cluster
          # j.name<-unlist(strsplit(j,"_mesh.ply"))
     
          rm.ply<-strsplit(k,".ply")
          rm.ply.underbar<-strsplit(unlist(rm.ply),"_")
          str.names<-(unlist(rm.ply.underbar))
          
          # insert names
          position.cluster<-which(visib_folder.select==k)
          R.result[position.cluster,1]<-str.names[3]
          R.result[position.cluster,2]<-str.names[1]
          R.result[position.cluster,3]<-str.names[2]
          
         # Calculate distance from the centroid to the edges
          BC<-vcgBary(oyster.cluster)
              BC.d<-as.data.frame(BC)
              X.BCm<-mean(BC.d$V1)
              Y.BCm<-mean(BC.d$V2)
              Z.BCm<-mean(BC.d$V3)
          centroid.distance<-max(sqrt((BC.d$V1-X.BCm)^2+(BC.d$V2-Y.BCm)^2+(BC.d$V3-Z.BCm)^2)) # Find the maximum distance from the center to the edge I in mm)
          SA<-vcgArea(oyster.cluster)
        
          R.result[position.cluster,4]<-centroid.distance
          R.result[position.cluster,5]<-SA
          
          # calculate roughness
              # import txt file
              setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist1/convert_name_full/")
          
              roughness.cluster = read.csv(paste0(unlist(rm.ply),".txt"),header=TRUE)
                 
              Rough<-(roughness.cluster$Roughness..4.01.);
              R.mean<-(mean(Rough,na.rm=TRUE))# 
              R.sd<-(sd(Rough,na.rm=TRUE))# 
              
              R.result[position.cluster,6]<-R.mean
              R.result[position.cluster,7]<-R.sd
              
      #return to original working directory
      setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist1")

      }   
      
#  Print out Summarise Table
R.result.d<-as.data.frame(R.result)
str(R.result.d)

# change to numerics
R.result.d$Cluster<-as.character(R.result.d$Cluster)
R.result.d$POV<-as.character(R.result.d$POV)
R.result.d$dist.sphere<-as.numeric(as.character(R.result.d$dist.sphere))
R.result.d$Centroid.distance<-as.numeric(as.character(R.result.d$Centroid.distance))
R.result.d$SA<-as.numeric(as.character(R.result.d$SA))
R.result.d$R.mean<-as.numeric(as.character(R.result.d$R.mean))
R.result.d$R.SD<-as.numeric(as.character(R.result.d$R.SD))

str(R.result.d)

# Data exploratory 
  R.result.d %>%
    group_by(POV, dist.sphere) %>%
    summarise(Surface =mean(SA), SurfaceSD=sd(SA), roughness=mean(R.mean), roughnessSD=sd(R.SD))
   

  R.result.d %>%
    group_by(Cluster) %>%
    summarise(Surface =mean(SA), SurfaceSD=sd(SA), roughness=mean(R.mean), roughnessSD=sd(R.SD))
  
# Centroid.distance
  Cent.dist<-R.result.d %>% 
    group_by(Cluster) %>%
    summarise(distance=max(Centroid.distance))
  mean(Cent.dist$distance); sd(Cent.dist$distance)
  summary(Cent.dist$distance)
  
  
# Cluster Calculate SA difference
  # Read surface and roughness of full CT images
  setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/")

  (fullcluster_folder <- list.files(pattern = ".ply"))
  fullcluster_folder.select<-fullcluster_folder # [c(41:45,81:86,146:150,211:216)]
  (fullcluster_folder.select.order <-mixedsort(sort(fullcluster_folder.select))) # use gtools to sort in number
  
  # make result saving matrix
  fullcluster.result<-matrix(0,length(fullcluster_folder.select.order),4)
  rownames(fullcluster.result)<-strsplit(fullcluster_folder.select.order,"_mesh.ply");fullcluster.result
  colnames(fullcluster.result)<-c("Cluster","SA","R.mean","R.SD")
  

        
          for (k in fullcluster_folder.select.order){ # "A160_3_Cluster39.ply"){ # 
            #read ply file 
            setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/")
            oyster.cluster<-vcgPlyRead(k)
            
            # full.cluster<-ply2mesh("Cluster4_mesh.ply") # only for # 4 cluster
            # j.name<-unlist(strsplit(j,"_mesh.ply"))
            
            rm.ply<-strsplit(k,"_mesh.ply")
            str.names<-(unlist(rm.ply))
            
            # insert names
            position.cluster<-which(fullcluster_folder.select.order==k)
            fullcluster.result[position.cluster,1]<-str.names
         
            
            # Surface area
            SA<-vcgArea(oyster.cluster)
            fullcluster.result[position.cluster,2]<-SA
            
            # calculate roughness
            # import txt file
            setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Roughness_fullcluster/")
            
            roughness.cluster = read.csv(paste0(unlist(rm.ply),".txt"),header=TRUE)
            
            Rough<-(roughness.cluster$Roughness..4.01.);
            R.mean<-(mean(Rough,na.rm=TRUE))# 
            R.sd<-(sd(Rough,na.rm=TRUE))# 
            
            fullcluster.result[position.cluster,3]<-R.mean
            fullcluster.result[position.cluster,4]<-R.sd
            
            #return to original working directory
            setwd("/Volumes/Seagate Backup Plus Drive/Full cluster_mesh/Visible_view160_dist1")
            
          }   
  
 
  fullcluster.result.d<-as.data.frame(fullcluster.result)
  # change to numerics
  fullcluster.result.d$Cluster<-as.character(fullcluster.result.d$Cluster)
  fullcluster.result.d$SA<-as.numeric(as.character(fullcluster.result.d$SA))
  fullcluster.result.d$R.mean<-as.numeric(as.character(fullcluster.result.d$R.mean))
  fullcluster.result.d$R.SD<-as.numeric(as.character(fullcluster.result.d$R.SD)) 
   
  fullcluster.result.d
   
  R.result.d
  
  
  
  ####################################################################################
  #### 3. Summarize results
  
  ## Summmarize the entire results
  entire.result<-matrix(0,54,11)
  rownames(entire.result)<-strsplit(fullcluster_folder.select.order,"_mesh.ply");entire.result
  colnames(entire.result)<-c("Cluster","SA_CT","SA_SFM1.160","SA_SFM1.80","SA_SFM3.160","SA_SFM3.80","RG_CT","RG_SFM1.160","RG_SFM1.80","RG_SFM3.160","RG_SFM3.80")
  
  for (m in fullcluster_folder.select.order){
    # m="Cluster2_mesh.ply"
    rm.full.ply<-strsplit(m,"_mesh.ply")
  
    # surface area
  Full.SA<-fullcluster.result.d$SA[fullcluster.result.d$Cluster == rm.full.ply]
  SA.160_1<-R.result.d$SA[ R.result.d$Cluster ==rm.full.ply & R.result.d$POV =="A160" & R.result.d$dist.sphere =="1"]
  SA.160_3<-R.result.d$SA[ R.result.d$Cluster ==rm.full.ply & R.result.d$POV =="A160" & R.result.d$dist.sphere =="3"]
  SA.80_1<-R.result.d$SA[ R.result.d$Cluster ==rm.full.ply & R.result.d$POV =="A80" & R.result.d$dist.sphere =="1"]
  SA.80_3<-R.result.d$SA[ R.result.d$Cluster ==rm.full.ply & R.result.d$POV =="A80" & R.result.d$dist.sphere =="3"]
  
  position.full.cluster<-which(fullcluster_folder.select.order==m)
  entire.result[position.full.cluster,1]<-unlist(rm.full.ply)
  
  entire.result[position.full.cluster,2]<-Full.SA/SA.160_1
  entire.result[position.full.cluster,3]<-1
  entire.result[position.full.cluster,4]<-SA.80_1/SA.160_1
  entire.result[position.full.cluster,5]<-SA.160_3/SA.160_1
  entire.result[position.full.cluster,6]<-SA.80_3/SA.160_1
  
  # roughness
  Full.RG<-fullcluster.result.d$R.mean[fullcluster.result.d$Cluster == rm.full.ply]
  RG.160_1<-R.result.d$R.mean[ R.result.d$Cluster ==rm.full.ply & R.result.d$POV =="A160" & R.result.d$dist.sphere =="1"]
  RG.160_3<-R.result.d$R.mean[ R.result.d$Cluster ==rm.full.ply & R.result.d$POV =="A160" & R.result.d$dist.sphere =="3"]
  RG.80_1<-R.result.d$R.mean[ R.result.d$Cluster ==rm.full.ply & R.result.d$POV =="A80" & R.result.d$dist.sphere =="1"]
  RG.80_3<-R.result.d$R.mean[ R.result.d$Cluster ==rm.full.ply & R.result.d$POV =="A80" & R.result.d$dist.sphere =="3"]
  
  entire.result[position.full.cluster,7]<-Full.RG
  entire.result[position.full.cluster,8]<-RG.160_1
  entire.result[position.full.cluster,9]<-RG.80_1
  entire.result[position.full.cluster,10]<-RG.160_3
  entire.result[position.full.cluster,11]<-RG.80_3
  
  
  }
  
  
  # change to numerics
  entire.result.d<-as.data.frame(entire.result)
  entire.result.d$Cluster<-as.character(entire.result.d$Cluster)
  entire.result.d$SA_CT<-as.numeric(as.character(entire.result.d$SA_CT))
  entire.result.d$SA_SFM1.160<-as.numeric(as.character(entire.result.d$SA_SFM1.160))
  entire.result.d$SA_SFM1.80<-as.numeric(as.character(entire.result.d$SA_SFM1.80))
  entire.result.d$SA_SFM3.160<-as.numeric(as.character(entire.result.d$SA_SFM3.160))
  entire.result.d$SA_SFM3.80<-as.numeric(as.character(entire.result.d$SA_SFM3.80))
  
  entire.result.d$RG_CT<-as.numeric(as.character(entire.result.d$RG_CT))
  entire.result.d$RG_SFM1.160<-as.numeric(as.character(entire.result.d$RG_SFM1.160))
  entire.result.d$RG_SFM1.80<-as.numeric(as.character(entire.result.d$RG_SFM1.80))
  entire.result.d$RG_SFM3.160<-as.numeric(as.character(entire.result.d$RG_SFM3.160))
  entire.result.d$RG_SFM3.80<-as.numeric(as.character(entire.result.d$RG_SFM3.80))
  
  entire.result.d
  
  # Final products
  fullcluster.result.d
  
  R.result.d
  
  entire.result.d

  #Create the Table for the paper
  str(entire.result.d)
  apply(entire.result.d[,c(2:11)],2,mean)
  apply(entire.result.d[,c(2:11)],2,sd)

  
  
  ####################################################################################
  #### 4. Compare using ANOVA and Posthoc
  
  
  #Compare roughness among methods
  
  library(tidyr)
  rough.entire<-entire.result.d[,7:11]
  rough.entire.gather<-gather(rough.entire,key=Methods, value=Roughness)
  res.aov<-aov(Roughness~Methods,data=rough.entire.gather)
  
  boxplot(rough.entire.gather$Roughness~rough.entire.gather$Methods)
  ggplot(rough.entire.gather,aes(x=Methods,y=Roughness))+
    geom_boxplot()+
    geom_jitter(shape=16, position=position_jitter(0.2))+
    scale_x_discrete(labels=c("RG_CT" = "CT", "RG_SFM1.160" = "HPR/1/160",
                              "RG_SFM1.80" = "HPR/1/80","RG_SFM3.160"="HPR/3/160","RG_SFM3.80"="HPR/3/80"))+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))
  
  summary(res.aov)
  TukeyHSD(res.aov)

  hist(rough.entire$RG_CT)
  hist(rough.entire$RG_SFM1.160)
  hist(rough.entire$RG_SFM1.80)
  hist(rough.entire$RG_SFM3.160)
  hist(rough.entire$RG_SFM3.80)
  
  plot(res.aov, 1)
  plot(res.aov, 2)
  
  #Compare Surface area among methods
  
  SA.entire<-entire.result.d[,2:6]
  SA.entire.gather<-gather(SA.entire,key=Methods, value=SurfaceArea)
  res.aov.SA<-aov(SurfaceArea~Methods,data=SA.entire.gather);summary(res.aov.SA)
 
  SA.entire.gather.outliarremove<-SA.entire.gather[SA.entire.gather$SurfaceArea<5,]
  
  ggplot(SA.entire.gather.outliarremove,aes(x=Methods,y=SurfaceArea))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(shape=16, position=position_jitter(0.2))+
    scale_x_discrete(labels=c("SA_CT" = "CT", "SA_SFM1.160" = "HPR/1/160",
                              "SA_SFM1.80" = "HPR/1/80","SA_SFM3.160"="HPR/3/160","SA_SFM3.80"="HPR/3/80"))+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold")) +
    ylab("Surface area ratio to HPR/1/160")
    
  
  unique()
  
  summary(res.aov.SA)
  TukeyHSD(res.aov.SA)
  
  hist(rough.entire$RG_CT)
  hist(rough.entire$RG_SFM1.160)
  hist(rough.entire$RG_SFM1.80)
  hist(rough.entire$RG_SFM3.160)
  hist(rough.entire$RG_SFM3.80)
  