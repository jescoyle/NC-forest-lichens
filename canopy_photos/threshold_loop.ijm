dir1 = getDirectory("Choose Source Directory "); 
list = getFileList(dir1); 
dir2 = getDirectory("Choose Directory to Save ");

thresholdTypes = newArray("Moments", "MaxEntropy", "IsoData")

for (j=0; j<thresholdTypes.length; j++)
{

t = thresholdTypes[j];
dir3 = dir2+t+File.separator;
File.makeDirectory(dir3);

for (i=0; i<list.length; i++) 
        { 
        
        if (File.isDirectory(dir1+list[i])){}	
        else{ 
                
                path = dir1+list[i]; 
                if (endsWith(path, ".db")){} 
                else{ 
                        
                        open(path); 
                        if (bitDepth!=24){}   
                        else	{ 
                                setBatchMode(true); 
                                title = File.nameWithoutExtension;
                                
			makeOval(384, 6, 1860, 1884);
			run("Fit Circle");
			run("Crop");
			run("Split Channels");
			run("Threshold...");
			setAutoThreshold(t);
			run("Convert to Mask");
	
			saveAs("Tiff", dir3+title+"_bw_"+t+".tif"); 
			close();
			close();
			close();
			setBatchMode(false); 

                                } 
                        } 
                } 
        }
}
