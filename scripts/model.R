library(here)
library(tidyverse)
library(patchwork)
library(pracma)
library(beepr)

addTaskCallback(function(...) {set.seed(42);TRUE})
options(warn=-1)
options(dplyr.summarise.inform = FALSE)

# Mee Simulation ----------------------------------------------------------

## Parameters:

NUM.reps <- 1 # The number of replicate simulations to run
## XX years total
NUM.gens.pre.reserve <- 12 # The number of generations before any fishery
NUM.gens.pre.fishing <- 6 # The number of generations of fishing before reserves are installed
NUM.gens.post.fishing <- 6 # The number of generations with the reserve installed
years = NUM.gens.pre.reserve+NUM.gens.pre.fishing+NUM.gens.post.fishing

NS.patches <- 10 # the number of patches on the north-south axis
EW.patches <- 10 # the number of patches on the east-west axis
patch.size <- 100 # the width and height of each grid cell in nautical miles (COULD BE METERS?)
## View the "world" coordinates:
view.world <- array(seq(1,NS.patches*EW.patches),c(NS.patches,EW.patches))
view.world

sb <- 0.37 # survival proportion for babies
s <- 0.37 # survival proportion
dd <- 0.0005 # density dependence of baby survival 
fecundity <- 2000 # The number of babies produced, on average, by each adult female each year.
maturity.age <- 4 # The average age at which individuals mature (i.e., the age at which 50% of individuals are mature)
fished <- 0.5
buffer.fished <- 0 #buffer fishing pressure (lower than total = buffer zone, higher than total = fishing the line)
reserves.at <- c(33,43,53,63,23,
                 34,44,54,64,24,
                 35,45,55,65,25,
                 36,46,56,66,26,
                 37,47,57,67,27) # This determines which patches are marine reserves. Should be a list: e.g., for one reserve, c(369,370,371,372,389,390,391,392,409,410,411,412,429,430,431,432)
buffer.at <- c()
mover.distance <- 200 # Individuals move this distance on average every year

############################################################################
## Create the world

NUM.age.classes <- 3 #babies, juvenile, adult
NUM.sexes <- 2 #female male

world <- array(0, c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes))

############################################################################
## This populates the world.

init <- function() {
  init <- 100
  pop <- world
  pop[,,,] <- init
  return(pop)
}

############################################################################
## This function creates an array to tell the simulation the reserve locations

where.reserves <- function(reserves.at) {
  reserve.patches <- array(0, c(NS.patches, EW.patches))
  for(i in 1:length(reserves.at)) {
    x <- ((reserves.at[i]-1) %/% NS.patches) + 1
    y <- ((reserves.at[i]-1) %% NS.patches) + 1
    reserve.patches[y,x] <- 1
  }
  return(reserve.patches)
}
reserve.patches <- where.reserves(reserves.at)

############################################################################
## This function creates an array to tell the simulation the buffer/fishing the line locations

where.buffer <- function(buffer.at) {
  buffer.patches <- array(0, c(NS.patches, EW.patches))
  for(i in 1:length(buffer.at)) {
    x <- ((buffer.at[i]-1) %/% NS.patches) + 1
    y <- ((buffer.at[i]-1) %% NS.patches) + 1
    buffer.patches[y,x] <- 1
  }
  return(buffer.patches)
}
buffer.patches <- where.buffer(buffer.at)

############################################################################
## This function initializes the habitat grid 

init.habitat = function() {
  
  
  
}


############################################################################
## This function causes adults to reproduce in spawning areas


spawn <- function(pop) {
  
  fec <- fecundity
  
  # All females produce the same mean number of eggs
  NUM.eggs <- Reshape(rpois(NS.patches * EW.patches,fec*pop[,,3,1]), NS.patches, EW.patches)
  # Divide zygotes 50:50 among the sexes
  babies.f <- rbinom(NS.patches * EW.patches,NUM.eggs,0.5)
  babies.m <- NUM.eggs - babies.f
  # Female babies
  pop[,,1,1] <- pop[,,1,1] + Reshape(babies.f, NS.patches, EW.patches)
  # Male babies
  pop[,,1,2] <- pop[,,1,2] + Reshape(babies.m, NS.patches, EW.patches)
  return(pop)
}

############################################################################
## This function determines natural survival and recruitment within each grid cell.

p <- 1/(maturity.age)

recruit <- function(pop) {
  recruit.array <- world
  # Some babies survive and recruit to juvenile age class
  s1 <- sb
  s <- s
  for(j in 1:NUM.sexes) {
    recruit.array[,,1+1,j] <- recruit.array[,,1+1,j] + Reshape(rbinom(NS.patches * EW.patches,pop[,,1,j],s1), NS.patches, EW.patches)
  }
  # Some juveniles survive
  for(j in 1:NUM.sexes) {
    juvies.surviving <- Reshape(rbinom(NS.patches * EW.patches,pop[,,2,j],s), NS.patches, EW.patches)
    # Some juveniles recruit to adult age class
    juvies.recruiting <- Reshape(rbinom(NS.patches * EW.patches,juvies.surviving,p), NS.patches, EW.patches)
    juvies.staying <- juvies.surviving-juvies.recruiting
    recruit.array[,,2+1,j] <- recruit.array[,,2+1,j] + juvies.recruiting
    # The rest of the juveniles remain in the juvenile age class
    recruit.array[,,2,j] <- recruit.array[,,2,j] + juvies.staying
  }
  # Some adults survive
  for(j in 1:NUM.sexes) {
    recruit.array[,,3,j] <- recruit.array[,,3,j] + Reshape(rbinom(NS.patches * EW.patches,pop[,,3,j],s), NS.patches, EW.patches)
  }
  return(recruit.array)
}

############################################################################
## This function determines fishing mortality within each grid cell, depending whether the cell is a reserve.

fishing <- function(pop,gen) {
  if(gen <= pre.reserve.gens+pre.fishing.gens) {
    each.patch.pop <- array(0,c(NS.patches,EW.patches))
    for(i in 2:NUM.age.classes) {
      for(j in 1:NUM.sexes) {
        each.patch.pop[,] <- each.patch.pop[,] + pop[,,i,j]
      }
    }
    mean.per.patch.pop <- mean(each.patch.pop)
    ff <- mean.per.patch.pop*(1/fished-1)
    patch.pop <- rowSums(pop[,,c(2,3),], dims = 2)
    f <- patch.pop/(ff+patch.pop)
    for(i in 2:NUM.age.classes) {
      for(j in 1:NUM.sexes) {
        pop[,,i,j] <- Reshape(rbinom(NS.patches * EW.patches,pop[,,i,j],(1-f)), NS.patches, EW.patches)
      }
    }
  }
  if(gen > pre.reserve.gens+pre.fishing.gens) {
    reserve.area <- sum(reserve.patches)/(NS.patches*EW.patches)
    buffer.area <- sum(buffer.patches)/(NS.patches*EW.patches)
    fished.adj <- (fished - (buffer.area*buffer.fished)) * 1/(1-(reserve.area + buffer.area))
    each.patch.pop <- array(0,c(NS.patches,EW.patches))
    for(i in 2:NUM.age.classes) {
      for(j in 1:NUM.sexes) {
        each.patch.pop[,] <- each.patch.pop[,] + pop[,,i,j]
      }
    }
    each.patch.pop = ifelse(reserve.patches == 1 | buffer.patches == 1, NaN, each.patch.pop)
    mean.per.patch.pop <- mean(each.patch.pop,na.rm=TRUE)
    ff <- mean.per.patch.pop*(1/fished.adj-1)
    patch.pop <- rowSums(pop[,,c(2,3),], dims=2)
    patch.pop = ifelse(reserve.patches == 1 | buffer.patches == 1, NaN, patch.pop)
    f <- patch.pop/(ff+patch.pop)
    if(buffer.fished != 0) {
      f = ifelse(buffer.patches == 1, buffer.fished, f)
    }
    for(i in 2:NUM.age.classes) {
      for(j in 1:NUM.sexes) {
        fished.array = Reshape(rbinom(NS.patches * EW.patches,pop[,,i,j],(1-f)), NS.patches, EW.patches)
        fished.array[is.na(fished.array)] = 0
        pop[,,i,j] <- ifelse(is.na(fished.array[,]),pop[,,i,j],fished.array[,])
      }
    }
  }
  return(pop)
}

############################################################################
## This function determines how far each individual moves. Movement distance for each genotype is drawn from a negative bimonial function. Babies do not move between grid cells. 

move <- function(pop) {
  
  movers <- mover.distance # Individuals with move this distance on average, in nautical miles

  move.array <- world
  
  for(lat in 1:NS.patches) {
    for(lon in 1:EW.patches) {
      for(i in 2:NUM.age.classes) {
        for(j in 1:NUM.sexes) {
          mean.dist = movers
          if(pop[lat,lon,i,j] > 0) {
            # movers are subtracted from the present grid cell
            move.array[lat,lon,i,j] <- move.array[lat,lon,i,j] - pop[lat,lon,i,j]
            # determine the distribution of movement distances in nautical miles:
            dist <- rnbinom(pop[lat,lon,i,j], mu = mean.dist, size = 1)
            # determine the direction of each move
            theta <- runif(pop[lat,lon,i,j],0,2*pi)
            # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
            # f.adj <- function(x, u) x-cos(x)*sin(x) - u
            # my.uniroot <- function(x) uniroot(f.adj, c(0, 2*pi), tol = 0.0001, u = x)$root
            # theta <- vapply(theta, my.uniroot, numeric(1))
            # convert direction and distance into a distance in the x-direction (longitude)
            x <- cos(theta)*dist
            # bounce off edges (assume fish start in centre of cell)
            for(m in 1:length(x)) {
              miss_x.edges <- FALSE
              while(miss_x.edges==FALSE) {
                if(x[m] <= patch.size*(EW.patches-lon)+patch.size/2) {
                  if(x[m] >= patch.size/2-lon*patch.size) {
                    miss_x.edges <- TRUE }
                }
                if(x[m] > patch.size*(EW.patches-lon)+patch.size/2) {
                  x[m] <- -(x[m]-2*(patch.size*(EW.patches-lon)+patch.size/2))
                  # distance penalty for hitting an edge
                  #x[m] <- x[m] + 1
                }
                if(x[m] < patch.size/2-lon*patch.size) {
                  x[m] <- -(x[m]-2*(patch.size/2-lon*patch.size))
                  # distance penalty for hitting an edge
                  #x[m] <- x[m] - 1
                }
              }
            }
            # convert direction and distance into a distance in the y-direction (latitude)
            y <- sin(theta)*dist
            # bounce off edges (assume fish start in centre of cell)
            for(m in 1:length(y)) {
              miss_y.edges <- FALSE
              while(miss_y.edges==FALSE) {
                if(y[m] <= patch.size*(NS.patches-lat)+patch.size/2) {
                  if(y[m] >= patch.size/2-lat*patch.size) {
                    miss_y.edges <- TRUE }
                }
                if(y[m] > patch.size*(NS.patches-lat)+patch.size/2) {
                  y[m] <- -(y[m]-2*(patch.size*(NS.patches-lat)+patch.size/2))
                  # distance penalty for hitting an edge
                  #y[m] <- y[m] + 1
                }
                if(y[m] < patch.size/2-lat*patch.size) {
                  y[m] <- -(y[m]-2*(patch.size/2-lat*patch.size))
                  # distance penalty for hitting an edge
                  #y[m] <- y[m] - 1
                }
              }
            }
            # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
            x.round = round(x/patch.size)
            y.round = round(y/patch.size)
            xy = as.data.frame(cbind(x.round,y.round))
            xy_sum = xy %>% 
              group_by(x.round, y.round) %>% 
              summarise(sum = n())
            xy_pivot = xy_sum %>% 
              ungroup() %>% 
              pivot_wider(values_from = sum, names_from = x.round)
            final_xy = xy_pivot %>% 
              select(-y.round)
            final_xy[is.na(final_xy)] <- 0
            final_xy = as.data.frame(final_xy)
            rownames(final_xy) = sort(unique(xy_sum$y.round))
            # populate the move.array with movers (and stayers)
            for(xx in 1:length(unique(xy_sum$x.round))) {
              for(yy in 1:length(unique(xy_sum$y.round))) {
                move.array[lat+as.numeric(row.names(final_xy)[yy]),lon+as.numeric(names(final_xy)[xx]),i,j] <- move.array[lat+as.numeric(row.names(final_xy)[yy]),lon+as.numeric(names(final_xy)[xx]),i,j] + final_xy[yy,xx]
              }
            }
          }
        }
      }
    }
  }
  # add move.array to pop to finish movement
  return(pop+move.array)
}

############################################################################
## THE SIMULATION ##########################################################

reps <- NUM.reps

pre.reserve.gens <- NUM.gens.pre.reserve
pre.fishing.gens <- NUM.gens.pre.fishing
post.fishing.gens <- NUM.gens.post.fishing
gens <- pre.reserve.gens+pre.fishing.gens+post.fishing.gens

output.array <- array(0 ,c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, gens, reps))

start_time <- Sys.time()

for(rep in 1:reps) {
  print(rep)
  pop <- init()
  for(t in 1:gens) {
    output.array[,,,,t,rep] <- pop
    pop <- spawn(pop)
    pop <- recruit(pop)
    if(t > pre.reserve.gens + pre.fishing.gens | t < pre.reserve.gens) {
      gen <- t
      pop <- fishing(pop,gen)
    }
    pop <- move(pop)
    print(pop)
    print(t)
  }
  gc() #clear memory
}
gc()

end_time <- Sys.time()
end_time - start_time

beepr::beep(5)

# Allie Explore -----------------------------------------------------------

# Output results into a dataframe
output_df = data.frame() #create dataframe to hold results

for(a in 1:reps) {
  world_sub <- array(0, c(NS.patches, EW.patches))
  for(b in 1:gens) {
    for(d in 1:NUM.sexes) {
      for(e in 1:NUM.age.classes) {
        world_sub = output.array[,,e,d,b,a] %>% 
          as.data.frame()
        world_sub$rep = paste0(a)
        world_sub$generation = paste0(b)
        world_sub$sex = paste0(d)
        world_sub$age = paste0(e)
        world_sub$lat = c(1:NS.patches)
        output_df = bind_rows(output_df, world_sub)
      }
    }
  }
}


# Wrangle dataframe into plottable format
output_df = output_df %>% 
  pivot_longer(V1:V10,
               names_to = "lon",
               values_to = "pop") %>% 
  mutate(lon = case_when(
    lon == "V1" ~ 1,
    lon == "V2" ~ 2,
    lon == "V3" ~ 3,
    lon == "V4" ~ 4,
    lon == "V5" ~ 5,
    lon == "V6" ~ 6,
    lon == "V7" ~ 7,
    lon == "V8" ~ 8,
    lon == "V9" ~ 9,
    lon == "V10" ~ 10
  )) %>%
  mutate(sex = case_when(
    sex == 1 ~ "female",
    sex == 2 ~ "male"
  )) %>% 
  mutate(sex = as.factor(sex)) %>%
  mutate(age = case_when(
    age == 1 ~ "baby",
    age == 2 ~ "juvenile",
    age == 3 ~ "adult"
  )) %>% 
  mutate(age = as.factor(age)) %>%
  mutate(lat = as.numeric(lat)) %>% 
  mutate(lon = as.numeric(lon))


#Summarize pop size and frequency by genotype
output_sum = output_df %>% 
  group_by(lat, lon, rep, generation) %>% 
  summarise(pop_sum = sum(pop)) 

#write_csv(output_sum, here("output", ".csv"))

# output_sum = read_csv(here("output", "3x3NoClimate8F.csv"))

plot_sum = output_sum %>% 
  filter(generation %in% c(12,14,16,18,20,22,24)) %>% 
  mutate(generation = as.numeric(generation))

plot = ggplot(plot_sum, aes(lon, lat, fill = pop_sum)) +
  geom_tile() + 
  facet_grid(~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size") +
  theme_bw() +
  scale_fill_gradient2(low = "gainsboro", high = "midnightblue", mid = "skyblue3", midpoint = 2)

plot

#ggsave(plot, file=paste0(".pdf"), path = here("figs"), height = 11, width = 8)

