
# ------------------------------------------------------------------------------------
# Functions for use in coexistence calcualtions

# Determine equilibrium conditions for each species in isolation 
pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0)
  return(N)
}

# invader population growth rate one time step forward
pop.invade <- function (N0, resident, s, g, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_inter*resident)
  return(N)
}

# resident population growth rate one time step forward
pop.resident <- function (N0, resident, s, g, a_intra, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*resident + resident*(lambda*g)/(1+a_intra*resident+a_inter*N0)
  return(N)
}

# ------------------------------------------------------------------------------------
# run models
# first determine resident equilibrium abundances and low density growth rates
# without partitioning coexistence

bromus <- subset(model.dat, species=="Bromus")
plantago <- subset(model.dat, species=="Plantago")

treatments <- unique(plerbrho$waterN_treatment)

## Set germination and survival fractions from the literature
ps <- .6 # lauren guess
pg <- .92 # gulmon
bs <- .2 # lauren guess
bg <- .98 # gulmon

N0 <- 100
time <- 120
# N_bromus <- rep(NA, time)
# N_bromus[1] <- N0


grwrdat <- data.frame(waterN_treatment = as.character(), bromus_epsilon_0 = as.numeric(), plantago_epsilon_0 = as.numeric(), resident_bromus_epsilon_0 = as.numeric(), 
                      resident_plantago_epsilon_0 = as.numeric(), plantago_grwrChesson = as.numeric(), bromus_grwrChesson = as.numeric(),
                      plantago_resident = as.numeric(), bromus_resident = as.numeric())

for (i in 1:length(treatments)){

  bromus_no_var <- rep(NA, time)
  bromus_no_var[1] <- N0
  
for (t in 1:time) {
  bromus_no_var[t+1] <- pop.equilibrium(N0=bromus_no_var[t], s=bs, g=bg, 
                                       a_intra=bromus$aiB[i], lambda=bromus$lambda[i])
}



# check output
plot(seq(1:(time+1)), bromus_no_var, type="l")


# for plantago
N0 <- 70
plantago_no_var <- rep(NA, time)
plantago_no_var[1] <- N0

for (t in 1:time) {
  plantago_no_var[t+1] <- pop.equilibrium(N0=plantago_no_var[t], s=ps, g=pg,  
                                          a_intra=plantago$aiP[i], lambda=plantago$lambda[i])
}

# check output
plot(seq(1:(time+1)), plantago_no_var, type="l")

# now invade each species into the resident
bromus_invade_no_var <- pop.invade(N0=1, resident=plantago_no_var[time], s=bs, g=bg, 
                                   a_inter=bromus$aiP[i], lambda=bromus$lambda[i])

plantago_invade_no_var <- pop.invade(N0=1, resident=plantago_no_var[time], s=ps, g=pg, 
                                     a_inter=plantago$aiB[i], lambda=plantago$lambda[i])

# determine any changes in the residents' abundances
bromus_resident_no_var_next <- pop.resident(N0=1, resident=bromus_no_var[time], s=bs, g=bg, 
                                           a_intra = bromus$aiB[i],
                                           a_inter=bromus$aiP[i], lambda=bromus$lambda[i])

bromus_resident_no_var <- bromus_no_var[time]/bromus_resident_no_var_next

# plantago
plantago_resident_no_var_next <- pop.resident(N0=1, resident=plantago_no_var[time], s=ps, g=pg, 
                                             a_intra = plantago$aiP[i],
                                             a_inter=plantago$aiB[i], lambda=plantago$lambda[i])

plantago_resident_no_var <- plantago_no_var[time]/plantago_resident_no_var_next

bromus_epsilon_0 <- log(bromus_invade_no_var)
plantago_epsilon_0 <- log(plantago_invade_no_var)
resident_bromus_epsilon_0 <- log(bromus_resident_no_var)
resident_plantago_epsilon_0 <- log(plantago_resident_no_var)
plantago_resident=plantago_no_var[time]
bromus_resident=bromus_no_var[time]


plantago_grwrChesson = plantago_epsilon_0-resident_bromus_epsilon_0
bromus_grwrChesson = bromus_epsilon_0-resident_plantago_epsilon_0
waterN_treatment = bromus$treatment[i]
tempout <- data.frame(waterN_treatment, bromus_epsilon_0, plantago_epsilon_0, 
                      resident_bromus_epsilon_0, resident_plantago_epsilon_0, plantago_grwrChesson, bromus_grwrChesson,
                      plantago_resident, bromus_resident)
grwrdat <- rbind(grwrdat, tempout)
}

grwrgraph <- grwrdat %>%
  select(waterN_treatment, plantago_grwrChesson, bromus_grwrChesson) %>%
  separate(waterN_treatment, c("water", "N"), sep = "_") %>%
  gather( "species", "grwrChesson", plantago_grwrChesson:bromus_grwrChesson) %>%
  separate(species, c("species", "todelete")) %>%
  select(-todelete)

ggplot(grwrgraph, aes(x = N, y = grwrChesson, fill = species)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~water)

ggplot(grwrgraph, aes(x = N, y = grwrChesson, fill = water)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~species)



residentgraph <- grwrdat %>%
  select(waterN_treatment, plantago_resident, bromus_resident) %>%
  separate(waterN_treatment, c("water", "N"), sep = "_") %>%
  gather( "species", "resident", plantago_resident:bromus_resident) %>%
  separate(species, c("species", "todelete")) %>%
  select(-todelete)

ggplot(residentgraph, aes(x = N, y = resident, fill = species)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~water)

ggplot(residentgraph, aes(x = N, y = resident, fill = water)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~species)
