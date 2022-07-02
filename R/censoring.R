de=function(x,y,w=1) 
{
  res= (x[1]-y[1])^2+(x[2]-y[2])^2
  if(length(x)==3)
  {
    res=res+w*w*(x[3]-y[3])^2
  }
  sqrt(res)
}

nearest<-function(p,points,w=1)
{
  res=1
  dmin=de(p,points[1,],w)
  if(nrow(points)>1)
  {
    for(i in 2:nrow(points))
    {
      di=de(p,points[i,],w)
      if(di<dmin)
      {
        res=i
        dmin=di
      }
    }
  }
  res 
}

createNA<-function(Data, m=2, p=0.75, schema=1,mesh_ref=NULL,locations=NULL,RDD_groups=12, mesh_ref_time=c(0,1),timelocations=NULL,w=1)
{
  if(schema==4) # ind in spazio, ind in tempo
  {
    data_vector=as.vector(Data)
    size=length(data_vector)
    data_vector[sample(1:size,size-round(size*p))]=NA
    Data=matrix(data = data_vector,nrow = nrow(Data),ncol = ncol(Data))
  }
  if(schema==1) # ind in spazio, dep in tempo
  {
    nodes=mesh_ref_time
    nobs=RDD_groups*p
    pts=timelocations
    for(s in 1:nrow(Data))
    {
      Itime = sample(1:length(mesh_ref_time),RDD_groups)
      seeds= nodes[Itime]
      Iobs=sample(1:RDD_groups,nobs)
      Obsmarker=rep(0,RDD_groups)
      Obsmarker[Iobs]=1
      for(k in 1:length(pts))
      {
        if(Obsmarker[nearest(c(pts[k],0),cbind(seeds,rep(0,length(seeds))))]==0)
        {
          Data[s,k]=NA
        }
      }
    }
  }
  if(schema==2) # dep in spazio, indep in tempo
  {
    nodes=mesh_ref$nodes
    nobs=RDD_groups*p
    pts=locations
    
    for(t in 1:ncol(Data))
    {
      Ispace = sample(1:nrow(nodes),RDD_groups)
      seeds= nodes[Ispace,]
      Iobs=sample(1:RDD_groups,nobs)
      Obsmarker=rep(0,RDD_groups)
      Obsmarker[Iobs]=1
      
      for(k in 1:nrow(pts))
      {
        if(Obsmarker[nearest(pts[k,],seeds)]==0)
        {
          Data[k,t]=NA
        }
      }
    }
  }
  if(schema==5)
  {
    Datavec=as.vector(Data)
    nodes=mesh_ref$nodes
    Ispace = sample(1:nrow(nodes),RDD_groups)
    #Itime =  sample(1:length(mesh_ref_time),RDD_groups)
    seedstime=runif(RDD_groups,head(mesh_ref_time,1),tail(mesh_ref_time,1))
    seeds= cbind(nodes[Ispace,1],nodes[Ispace,2],seedstime)#mesh_ref_time[Itime])
    nobs=RDD_groups*p
    Obsmarker=rep(0,RDD_groups)
    Iobs=sample(1:RDD_groups,nobs)
    Obsmarker[Iobs]=1
    pts=cbind(rep(locations[,1],length(timelocations)),rep(locations[,2],length(timelocations)),rep(timelocations,each=nrow(locations)))
    for(k in 1:nrow(pts))
    {
      if(Obsmarker[nearest(pts[k,],seeds,w)]==0)
      {
        Datavec[k]=NA
      }
    }
    Data=matrix(data = Datavec, nrow = nrow(locations),ncol=length(timelocations))
    # plot3d(seeds[Iobs,],col="green",size=10)
    # points3d(seeds[-Iobs,],col="red",size=10)
    # points3d(pts[which(!is.na(as.vector(Data))),],col="grey",size=3)
    
  }
 Data
}
