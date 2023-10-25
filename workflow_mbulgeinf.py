import subprocess 
import os

for i in range(513):
    slice1 = i*1000
    slice2 = (i+1)*1000
    if i == 512:
        slice1 = i*1000
        slice2 = 512821
    
    command = "/scratch/alpine/mahu8801/blackhole_data/scripts/commands/%i_%i.sh"%(slice1,slice2)
    outfile = '/scratch/alpine/mahu8801/blackhole_data/data/%i_%i.dat'%(slice1,slice2)
    
    
    with open('/scratch/alpine/mahu8801/blackhole_data/scripts/commands/mbulgeinf.sh', 'r') as file :
        filedata = file.read()
        
    filedata = filedata.replace('--slice1=s1', '--slice1=%i'%(slice1))
    filedata = filedata.replace('--slice2=s2', '--slice2=%i'%(slice2))
    filedata = filedata.replace('s1_s2.dat', '%i_%i.dat'%(slice1,slice2))
        
    with open(command, 'w') as file:
        file.write(filedata)
            
    subprocess.run(["sbatch", command],capture_output=True)
