require(igraph)

bfs <- function(G,s,t){
  
  n = nrow(G)
  parent = integer(n)
  v = integer(n)
  v[s]=1
  q=c(s)
  
  while (length(q)){
    #print (q)
    u=q[1]
    q=q[-1]
    for (i in 1:n){
      if (v[i] == 0 && G[u,i] > 0){
        q=c(q,i)
        v[i]=1
        parent[i]=u
      }
    }
  }

  return(c(v[t],parent))

}

SS <- function(G,source,sink,sinkweight = NULL){
  n = nrow(G)
  m = ncol(G)
  plt = graph.adjacency(G,weighted = TRUE)
  plot(plt)
  
  MAXFLOW = 0
  PATHS = list()
  
  parent = integer(n)-1
  bf = bfs(G,source,sink)
  pathexist = bf[1]
  parent = bf[-1]
  
  if (is.null(sinkweight)){
    while (pathexist){
      s = sink
      path = Inf
      gnodes = c()
      
      while (s != source){
        gnodes = c(gnodes,s)
        path = min(path,G[parent[s],s])
        s = parent[s]
        
      }
      
      MAXFLOW = MAXFLOW + path
      gnodes = c(gnodes,source)
      gnodes = c(gnodes , path)
      
      PATHS = append(PATHS,list(gnodes))

      back = sink
      while (back != source)
      {
        temp = parent[back]
        G[temp,back] = G[temp,back] - path
        G[back,temp] = G[back,temp] + path
        back = parent[back]
      }
      bf = bfs(G,source,sink)
      pathexist = bf[1]
      parent = bf[-1]
    }
  }
  else{
    while (pathexist){
      s = sink
      path = Inf
      gnodes = c(source)
      
      while (s != source){
        
        gnodes = c(gnodes,s)
        path = min(path,G[parent[s],s],sinkweight)
        s = parent[s]
        
      }
      
      MAXFLOW = MAXFLOW + path
      gnodes = c(gnodes , path)
      PATHS = append(PATHS,list(gnodes))

      
      back = sink
      while (back!= source)
      {
        temp = parent[back]
        G[temp,back] = G[temp,back] - path
        G[back,temp] = G[back,temp] + path
        back = parent[back]
      }
      
      if (path == sinkweight) {
        break
      }
      else
      {
        sinkweight = sinkweight - path
      }
      
      bf = bfs(G,source,sink)
      pathexist = bf[1]
      parent = bf[-1]
    }
  }
  return (list(MAXFLOW,G,PATHS))

}

MaxFlow <- function (G,sources,sinks,sinkweights = NULL){
  GRAPH = G
  GRAPH = graph.adjacency(GRAPH , weighted = TRUE)
  nsource = length(sources)
  nsinks = length(sinks)
  
  if (!(is.null(sinkweights))){
    nweights = length(sinkweights)
    if (nsinks != nweights){
      print ("Number of sinks and sinkweights does not match ...")
      return (NULL)
    }
  }
  
  if (is.null(sinkweights)){
    if (nsinks == 1){
      
      MAXFLOWSFROM = integer(nsource)
      MAXFLOWSTO = 0
      GRAPHPATHS = list()
      
      for (i in 1:nsource)
      {
        ReturnValues = SS(G,sources[i],sinks)
        MAXFLOWSFROM[i] = ReturnValues[[1]]
        MAXFLOWSTO = MAXFLOWSTO + ReturnValues[[1]]
        G = ReturnValues[[2]]
        GRAPHPATHS = append(GRAPHPATHS,ReturnValues[[3]])
      }
      UniqueNodes = c()
      for (i in GRAPHPATHS){
        len = length(i)
        UniqueNodes = unique(c(UniqueNodes , i[-len]))
      }
      print(UniqueNodes)
      print (GRAPHPATHS)
      plot.igraph(GRAPH,mark.groups = UniqueNodes)
      cat("The Maximum flow is ",MAXFLOWSTO)
      return (GRAPHPATHS)
    }
    else{
      
      MAXFLOWS = matrix(0,nsource,nsinks)
      GRAPHPATHS = list()
      
      #for (i in 1:nsource){
      #  GRAPHPATHS = append(GRAPHPATHS , list(0))
      #}
      
      print (GRAPHPATHS)
      for(i in 1:nsource){
        tempPath = list()
        for (j in 1:nsinks){
          ReturnValues = SS(G,sources[i],sinks[j])
          MAXFLOWS[i,j] = ReturnValues[[1]]
          G = ReturnValues[[2]]
          tempPath = append(tempPath , ReturnValues[[3]])
        }
        GRAPHPATHS = append(GRAPHPATHS , tempPath)
      }
      UniqueNodes = c()
      for (i in GRAPHPATHS){
        len = length(i)
        UniqueNodes = unique(c(UniqueNodes , i[-len]))
      }
      print(UniqueNodes)
      
      plot.igraph(GRAPH,mark.groups = UniqueNodes)
      cat("The Maximum flow is ",MAXFLOWS)
      print (GRAPHPATHS)
      return (GRAPHPATHS)
      
    }
  }
  
  else{
    MAXFLOWS = matrix(0,nsource,nsinks)
    GRAPHPATHS = list()
    flag = TRUE
    checkflows = matrix(0,nsource,nsinks)
    
    while (TRUE){
      flag = TRUE
      for(j in 1:nsinks){
        tempPath = c()
        for(i in 1:nsource){
          ReturnValues = SS(G,sources[i],sinks[j],sinkweights[j])
          MAXFLOWS[i,j] = MAXFLOWS[i,j] + ReturnValues[[1]]
          
          if (sum(MAXFLOWS[,j]) >= sinkweights[j]){
            sinkweights[j] = Inf
          }
          
          if (ReturnValues[[1]] != 0){
            flag =FALSE
          } 
          checkflows[i,j] = ReturnValues[[1]]
          
          G = ReturnValues[[2]]
          tempPath = append(tempPath , ReturnValues[[3]])
        }
        GRAPHPATHS = append(GRAPHPATHS , tempPath)
      }
      if (flag) break
    }
    print (GRAPHPATHS)
    print (MAXFLOWS)
    
    UniqueNodes = c()
    for (i in GRAPHPATHS){
      len = length(i)
      UniqueNodes = unique(c(UniqueNodes , i[-len]))
    }
    print(UniqueNodes)
    
    plot.igraph(GRAPH,mark.groups = UniqueNodes)
    cat("The Maximum flow is ",MAXFLOWS)
    
    return (GRAPHPATHS)
  }
}

#G = matrix(c(0,16,13,0,0,0,0,0,10,12,0,0,0,4,0,0,14,0,0,0,9,0,0,20,0,0,0,7,0,4,0,0,0,0,0,0),6,byrow=TRUE)

a=MaxFlow(G,c(1,2),c(10,5),c(21,15))
