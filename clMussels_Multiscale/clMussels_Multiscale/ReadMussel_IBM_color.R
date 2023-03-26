# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Plotting script for the CUDA implementation of the                          #
# Multi-scale mussel bed patterns model of Liu et al 2014                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

remove(list=ls()) # Remove all variables from memory

on=1;off=0;

setwd("/Simulations/Opencl/clMussels_MultiScale/clMussels_MultiScale/Data")

require(fields)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Program settings and parameters
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Movie=on
Wait=off
SmallMovie=off
StraightToEnd=off

DPI=144

#Width = 400

Location=c(400,400)
Size = 400

WinWidth = 1440
WinHeight = 720
Resolution="1440x720"

# Graphical parameters & palette definitions
Algae.palette = colorRampPalette(c("Blue", "Green"))
mud.palette = colorRampPalette(c("aliceblue", "cornsilk", "cornsilk4", "darkgrey"))
mussel.palette = colorRampPalette(c("aliceblue", "cornsilk", "cornsilk4", "black"));

# The maximal value of P, O and W to be plotted. If one goes over that, the value is capped
AGraphMax = 2  
MGraphMax = 100
SGraphMax = 40

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Displaying the 2 plots
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

PlotBed = function (X,Y,popM,popS){

  par(mar=c(2.5, 1, 1.5, 1), mfrow=c(1,2))
  
  Dim=dim(popS)
  
  numInd=length(X)
  
  x=(1:Dim[1]-0.5)/Dim[1]*DomainSize
  y=(1:Dim[2]-0.5)/Dim[2]*DomainSize
  
  # First figure, setting margin, plotting, and adding a title
  LX=c(Location[1],Location[1],Location[1]+Size,Location[1]+Size,Location[1])
  LY=c(Location[2],Location[2]+Size,Location[2]+Size,Location[2],Location[2])
  
  plot(LX,LY, type='l', bty="n", cex=1, asp=1, xaxt='n', yaxt="n", col="red")
  
  clip(Location[1],Location[1]+Size,Location[2],Location[2]+Size)
  
  image(x=x, y=y, z=pmin(t(popS),SGraphMax), zlim=c(0,SGraphMax), yaxt="n", 
        xlim=c(Location[1],Location[1]+Size), ylim=c(Location[1],Location[1]+Size),
        add=TRUE, col = mud.palette(255),asp=1, bty="n", useRaster=T)
  
  image.plot(z=pmin(popS,SGraphMax), legend.only=T,
             zlim=c(0,SGraphMax), cex.axis=0.8,
             yaxt="n", horizontal = F, add=TRUE,
             col = mussel.palette(255),asp=1, bty="n", useRaster=T,
             legend.shrink = 0.93, legend.width = 1,smallplot= c(.93,.96,0.13,0.91))  
  
  rect(Location[1],Location[2],Location[1]+Size,Location[2]+Size, border="red")
  
  InSel = (X>=Location[1]&X<=Location[1]+Size&Y>=Location[2]&Y<=Location[2]+Size)
  Xs=X[which(InSel)];Ys=Y[which(InSel)];
  points(Xs,Ys, pch=20, cex=100/Size, bg=86)
  
  NumMussels=sprintf('Number within frame : %1.0f',length(Xs))
  title(main= NumMussels, line = -0.4, cex.main=0.9, font.main=1) 
  
  axis(side=1, line=-0.5, cex.axis=1, tck = -0.02, mgp=c(3, .5, 0))
  
  title(xlab="Colors indicate accumulated sediment (mm)", line=1, cex.lab=0.9)
  
  mtext(text=paste("Time : ",sprintf("%1.0f",(jj+1)/NumFrames*EndTime/24),
                   "of" ,sprintf("%1.0f",EndTime/24), "days"), 
        side=3, line=-1, cex=0.9, outer=TRUE)
  
  # Second figure, seting margin, plotting, adding title and counter text      
  LX=c(0,0,DomainSize,DomainSize,0)
  LY=c(0,DomainSize,DomainSize,0,0)
  
  plot(LX,LY, type='l', bty="n", xlab="n", ylab="", asp=1, 
       xlim=c(0,DomainSize), ylim=c(0,DomainSize),
       xaxt='n', yaxt="n",col='black') 
  
  #rect(0,0,DomainSize,DomainSize, border="red")
  
  if (numInd<500000){
    image(x=x, y=y, z=pmin(t(popS),SGraphMax), zlim=c(0,SGraphMax), yaxt="n", 
          add=TRUE, col = mud.palette(255),asp=1, bty="n", useRaster=T)  
    points(X,Y, pch=16, cex=0.1, bty="n", xlab="n", ylab="", asp=1, 
         xaxt='n', yaxt="n",col='black') 
  } else {
    image(x=x, y=y, z=pmin(t(popM),MGraphMax), zlim=c(0,MGraphMax), yaxt="n", 
          add=TRUE, col = mussel.palette(255),asp=1, bty="n", useRaster=T)  
  }
  
  rect(Location[1], Location[2],Location[1]+Size, Location[2]+Size, border="red")
  
  axis(side=1, line=-0.5, cex.axis=1, tck = -0.02, mgp=c(3, .5, 0))
  
  NumMussels=sprintf('Total number of mussels : %1.0f',N_in_Time[jj+1])
  title(main= NumMussels, line = -0.4, cex.main=0.9, font.main=1) 
  
  title(xlab="Space (cm)", line=1, cex.lab=0.9)
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Reading the parameters from the data file and declaring variables
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

NumFrames = length(dir())/2

N_in_Time = 1:NumFrames

if (StraightToEnd==on) StartNumber=NumFrames-1 else StartNumber=0

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Opening a window and starting the display loop
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if (Movie==off) 
  if(Sys.info()["sysname"]=="Darwin"){
    #quartz(width=WinWidth/DPI, height=WinHeight/DPI, dpi=DPI)
    X11(width=WinWidth/DPI*1.5, height=WinHeight/DPI*1.5, type="cairo", antialias="none", pointsize=12*1.5)
  } else {
    X11(width=WinWidth/DPI, height=WinHeight/DPI)
    quartz(width=WinWidth/DPI, height=WinHeight/DPI)
  }

for (jj in StartNumber:(NumFrames-1)){  # Here the time loop starts 
  
   # If a movie is to be made, frames are written as jpegs
   if (Movie==on)
      jpeg(filename = sprintf("../Images/Rplot%03d.jpeg",jj),
           width = WinWidth, height = WinHeight, 
           units = "px", pointsize = 24,
           quality = 100,
           bg = "white", res = NA,
           type = "quartz")  
   
   # Reading the data from the files
   IBMname = sprintf('IBM%i.dat',jj)
   PDEname = sprintf('PDE%i.dat',jj)
  
   FID = file(IBMname, "rb")
   NumMussels = readBin(FID, integer(), n = 1, endian = "little");
   DomainSize = readBin(FID, integer(), n = 1, endian = "little");
   EndTime = readBin(FID, integer(), n = 1, endian = "little");
   X=matrix(nrow=NumMussels, ncol=1, readBin(FID, numeric(), size=4, n = NumMussels, endian = "little"));
   Y=matrix(nrow=NumMussels, ncol=1, readBin(FID, numeric(), size=4, n = NumMussels, endian = "little"));
   close(FID)
   
   FID = file(PDEname, "rb")
   WPOP = readBin(FID, numeric(), size=4, n = 1, endian = "little");
   HPOP = readBin(FID, numeric(), size=4, n = 1, endian = "little");
   
   popA=matrix(nrow=HPOP, ncol=WPOP, byrow=T, readBin(FID, numeric(), size=4, n = WPOP*HPOP, endian = "little"));
   popM=matrix(nrow=HPOP, ncol=WPOP, byrow=T, readBin(FID, numeric(), size=4, n = WPOP*HPOP, endian = "little"));
   popS=matrix(nrow=HPOP, ncol=WPOP, byrow=T, readBin(FID, numeric(), size=4, n = WPOP*HPOP, endian = "little"));
   close(FID)
   
   # recording means values per recorded frame
   N_in_Time[jj+1] = length(X)
   
   PlotBed(X,Y,popM,popS)
   
   # Finishing JPEG, or updating graph
   if (Movie==on) dev.off() else { 
     dev.flush()
     dev.hold()
   }
   
   # For debugging, this lets you go frame by frame
   if (Wait==on){
     cat ("Press [enter] to continue, [q] to quit")
     line <- readline()
     if (line=='q'){ stop() }
   }   
     
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Building the movie via ffmpeg
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if (Movie==on) { 
  
   # Building command line to run ffmpeg
   InFiles=paste(getwd(),"/../Images/Rplot%03d.jpeg", sep="")
   OutFile="../Mussel3D_color.mp4"
  
   CmdLine=sprintf("ffmpeg -y -r 10 -i %s -s %s -c:v libx264 -pix_fmt yuv420p -b:v 10000k %s", 
                   InFiles, Resolution, OutFile)
   
   # Executing the command
   cmd = system(CmdLine)
  
   # Unhash to immediately display the movie
   # if (cmd==0) try(system(paste("open ", paste(getwd(),"Mussels_PDE.mp4"))))
} 

OnGoing=TRUE; dev.flush()

while(OnGoing){

  par(mar=c(2.5, 1, 1.5, 1), mfrow=c(1,2))
  
  options(locatorBell = FALSE)
  loc=locator(n=1)
  Location=c(loc$x-0.5*Size,loc$y-0.5*Size)
  
  print(paste(Location[1],Location[2]))
  
  PlotBed(X,Y,popM,popS)
  
  OnGoing = !is.null(loc)
}

system('say All ready')
