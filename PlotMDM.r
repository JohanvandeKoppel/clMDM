# First setup of the model
remove(list=ls()) # Remove all variables from memory

on=1;off=0;
setwd('/Simulations/OpenCL/clMDM/clMDM')

require(fields)

FID = file("MDM_Output.dat", "rb")

#---------- Settings ---------------------------------------------------
Movie=off
Wait=off
WinWidth = 960
WinHeight = 720
DPI=144

NX = readBin(FID, integer(), n = 1, endian = "little");
NY = readBin(FID, integer(), n = 1, endian = "little");
NumFrames = readBin(FID, integer(), n = 1, endian = "little")-1;
EndTime=readBin(FID, integer(), n = 1, endian = "little");

if (Movie==off)   
  # Pick the one you prefer below. X11 is fast, quartz looks nicer. 
  # Windows users should use "windows"
  #quartz(width=WinWidth/DPI, height=WinHeight/DPI, dpi=DPI)
  X11(width=WinWidth/DPI*1.5, height=WinHeight/DPI*1.5, 
      type="cairo", antialias="none", pointsize=12*1.5)

for (jj in 0:(NumFrames)){  # Here the time loop starts 
  
   if (Movie==on)
     tiff(filename = sprintf("Images/Rplot%03d.tiff",jj),
          width = WinWidth, height = WinHeight, 
          units = "px", pointsize = 24,
          compression="none",
          bg = "white", res = NA,
          type = "quartz")  
   
   Data_C = matrix(nrow=NY, ncol=NX, readBin(FID, numeric(), size=4, n = NX*NY, endian = "little"));
   
   par(mar=c(3, 4, 3, 4) + 0.1)
   
   image(Data_C, zlim=c(0,2), xaxt="n", yaxt="n",
         col = c("red", "black","green"),
         asp=1, bty="n", useRaster=TRUE)  
   
   title("Mussel Disturbance Model")   
   
   mtext(text='Green = mussels;', col = "darkgreen", side=1, line=0.3, adj=0.2, cex=0.9)
   mtext(text='Red = disturbed;', col = "red", side=1, line=0.3, adj=0.53, cex=0.9)
   mtext(text='Black = bare', col = "black", side=1, line=0.3, adj=0.8, cex=0.9)
   
   mtext(text=paste("Time : ",sprintf("%1d",floor(jj/NumFrames*EndTime)),
                    "of" ,sprintf("%1d",EndTime, "Timesteps")), 
                    side=1, adj=0.5, line=1.6, cex=0.9)
   
   if (Movie==on) dev.off() else { 
     dev.flush()
     Sys.sleep(0.0)
     dev.hold()
   }
   if (Wait==on){
     cat ("Press [enter] to continue, [q] to quit")
     line <- readline()
     if (line=='q'){ stop() }
   } 
}

close(FID)

if (Movie==on) { 
  
   InFiles=paste(getwd(),"/Images/Rplot%03d.tiff", sep="")
   OutFile="Mussels_MDM.mp4"
  
   print(paste(" building :", OutFile))
  
   CmdLine=sprintf("ffmpeg -y -r 25 -i %s -c:v libx264 -pix_fmt yuv420p -b:v 5000k %s", InFiles, OutFile)
   cmd = system(CmdLine)
  
   # if (cmd==0) try(system(paste("open ", paste(getwd(),"Mussels_MDM.mp4"))))
} 

system('say All ready')
