# Mussel Disturbance Model
# Johan van de Koppel - Januari 2005, update September 2013

# Following Guichard et at, 2003. Mussel Disturbance Dynamics: Signatures of Oceanographic
# Forcing from Local Interactions. Am. Nat. 2003. Vol. 161, pp. 889-904.

# First setup of the model
remove(list=ls()) # Remove all variables from memory

require("simecol")

# ----- Model parameters -----------------------------------------------------
Alpha0  = 0.8         # Disturbance probability if one of neighbours is disturbed
Alpha2  = 0.5         # Recovery probability coefficient
Lambda0 = 0.0002      # Probability that an occupied cell becoming disturbed by chance

# Simulation parameters
n       = 100         # Dimensions of the simulated landscape
EndTime = 500         # Number of timesteps in the whole simulation
NoFrames= 500         # Number of frames displayed during the entire simulation

WinWidth=5            #      - Width of the simulation window 
WinHeight=5           #      - Height of the simulation window

# Initialisation
wdist8=matrix(nrow=3,c(1,1,1,1,0,1,1,1,1))
DistNeighbourhood=wdist8      # Matrix of influence for disturbance
GrowthNeighbourhood=wdist8    # Matrix of influence for growth

Frac=0.02                     # Fraction of cells occupied
DFrac=0.0                     # Fraction of cells disturbed 

# Initial values of the Cell matrix
R=runif(n*n)                  # random number for each cell
Cells = matrix(ncol=n,data=1+(R<Frac)-2*(R<DFrac))

ii=EndTime                    # This is the plotcounter

#------- Setting up the figure ------------------------------------------
## Open a graphics window (Darwin stands for a Mac computer)
if (Sys.info()["sysname"]=="Darwin"){
  quartz(width=WinWidth, height=WinHeight, 
         title="Mussel disturbance model")
} else
  windows(width = WinWidth, height = WinHeight,
          title="Mussel disturbance model")

par(mfrow=c(1,1))

# ------------ The simulation loop ---------------------------------------
print(system.time(
for (Time in 1:EndTime){

    # For each cell, the number of disturbed neighbours (e.g., have value 0).
    # I add one to avoid miscalculation at the edge
    SumDist=matrix(nrow=n,neighbours(x=Cells+1,state=1,wdist=DistNeighbourhood))
    
    # For each cell, the number of occupied neighbours (e.g., have value 2)
    SumOccu=matrix(nrow=n,neighbours(x=Cells,state=2,wdist=GrowthNeighbourhood))
   
    Nu0=SumDist>0    # Nu0=1 if at least one of the neighbour cells in disturbed  
    N=SumOccu/8      # N is the number of neighbours as a fraction  

    R=matrix(nrow=n,data=runif(n*n)) # random number for each cell
    
    # If the random number is bigger than N*Alpha2, then there is colonisation 
    # (1->2, so one is added)
    Colonisation = (Cells==1)&(R<=(N*Alpha2))
    
    # If the cell is occupied (2), it can erode if:
    # Change = 2 -> 0, so 2 is substracted
    Erosion = ((Cells==2)&(R<=(Alpha0*Nu0+Lambda0)))*-2
    
    # If disturbed (0), it will become bare (1), and so 1 will be added
    BecomeBare=(Cells==0)
    
    # Combining it all the calculate the new cell state
    Cells = Cells + Colonisation + Erosion + BecomeBare

    # Graphic representation of the model every now and then
    if (ii>=EndTime/NoFrames){
       image(Cells, zlim=c(0,2), xaxt="n", yaxt="n",
             col = c("red", "black","green"),
             sub = "Green = mussels; Red = eroding; Black=bare",
             main="Spatial view",
             xlab=paste("Time : ",sprintf("%4.0f",Time),
                       "of" ,sprintf("%4.0f",EndTime), "days"), cex=1) 
       
       # The following two lines prevent flickering of the screen
       dev.flush() # Force the model to update the graphs
       dev.hold()  # Put all updating on hold       

       ii=0    # Resetting the plot counter 

    }
    
    ii=ii+1 # The plot counter is updated
    
} ))

   


