import subprocess 
import os

for i in range(47):
    slice1 = i*1000
    else:
        slice2 = (i+1)*1000
    
command = "/scratch/alpine/mahu8801/blackhole_data/scripts/commands/%i_%i.sh"%(slice1,slice2)
outfile = '/scratch/alpine/mahu8801/blackhole_data/data/%i_%i.dat'%(slice1,slice2)

if os.path.exists(outfile) == False:
    with open('/scratch/alpine/mahu8801/blackhole_data/scripts/commands/command.sh', 'r') as file :
        filedata = file.read()
        
    filedata = filedata.replace('--slice1=0', '--slice1=%i'%(slice1))
    filedata = filedata.replace('--slice2=1000', '--slice2=%i'%(slice2))
    filedata = filedata.replace('data/1000draws.dat', 'data/new/%i_%i.dat'%(slice1,slice2))
        
    with open(command, 'w') as file:
        file.write(filedata)
            
    subprocess.run(["sbatch", command],capture_output=True)
