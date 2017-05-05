track.shift <- function(shift="sp",tree, node){ #track shift="sp" or "ext"
  if (!is.numeric(node)){
    nodeold <- node
    node <- as.numeric(gsub("[^0-9]","",nodeold))
    print(paste0("we identified that '",nodeold, "' refers to node number ", node ))
  }
  edge <- tree$edge
  shiftspm <- tree$shiftsp
  shiftextm <- tree$shiftext
  
  ###################################################### inner functions ######################################
  # tracing back the shifts on speciation for one species until origin, name of the attributes tell the node which change happened
  trajectoryshiftsp  <- function (trace){
    #trace should indicate the edge number that is to be followed until the origin
    trajectory  <- NULL
    while ( length( which(edge[,2] == trace)) ){
      atual  <- which(edge[,2] == trace)
      add <- shiftspm[shiftspm[,"node"]==edge[atual,2],"strength"]
      attributes(add)$names <- shiftspm[shiftspm[,"node"]==edge[atual,2],"node"]
      trajectory  <- c(add, trajectory)
      trace  <- edge[atual,1]
    }
    return(trajectory)
  }
  # tracing back the shifts on extinction for one species until origin, name of the attributes tell the node which change happened
  trajectoryshiftext  <- function (trace){
    #trace should indicate the edge number that is to be followed until the origin
    trajectory  <- NULL
    while ( length( which(edge[,2] == trace)) ){
      atual  <- which(edge[,2] == trace)
      add <- shiftextm[shiftextm[,"node"]==edge[atual,2],"strength"]
      attributes(add)$names <- shiftextm[shiftextm[,"node"]==edge[atual,2],"node"]
      trajectory  <- c(add, trajectory)
      trace  <- edge[atual,1]
    }
    return(trajectory)
  }
  ############################################ end inner functions ######################################
  
  if (shift=="sp"){
    print("sp")
    trajectory <- trajectoryshiftsp(trace=node)
  }else if (shift=="ext"){
    print("ext")
    trajectory <- trajectoryshiftext(trace=node)
  }else{
    print("please choose either `sp` (tracking shifts on speciation) or `ext` (tracking shifts on extinction) for the shift parameter")
  }
  
  return(trajectory)
  
}
