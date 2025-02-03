
# – Calculate the dissimilarity to present natural cells within the same biogeographical realm####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

median.distances.from.PN <- function(scenario, realm, data, distance.matrix) {
  
  #Select scenarios
  data.current = data[which(data$Scenario == scenario),]
  data.pn = data[which(data$Scenario == "Present-Natural"),]
  
  # Select cells
  cells.current.selection = which(data.current$Biogeo.realm == realm)
  cells.pn.selection = which(data.pn$Biogeo.realm == realm)
  
  # Calculate median distances
  apply(distance.matrix[cells.current.selection,cells.pn.selection],1,median.without.na)
}


# – Calculate the dissimilarity to present natural cells in other biogeographical realms####
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

median.distance.to.other.realms <- function(scenario, realm, data, distance.matrix) {
  
  #Select scenarios
  data.temp = data[which(data$Scenario == scenario),]
  
  # Select cells
  cells.within.realm.selection = which(data.temp$Biogeo.realm == realm)
  cells.outside.realm.selection = which(!(data.temp$Biogeo.realm == realm))
  
  # Calculate median distances
  apply(distance.matrix[cells.within.realm.selection,cells.outside.realm.selection],1,median.without.na)
}
