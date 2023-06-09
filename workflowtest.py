import subprocess 

slice1 = 538000
slice2 = 538818
with open('/scratch/alpine/mahu8801/blackhole_data/scripts/commands/command.sh', 'r') as file :
    filedata = file.read()
        
filedata = filedata.replace('--slice1=0', '--slice1=%i'%(slice1))
filedata = filedata.replace('--slice2=1000', '--slice2=%i'%(slice2))
filedata = filedata.replace('1000draws', '%i_%i'%(slice1,slice2))
    
with open('/scratch/alpine/mahu8801/blackhole_data/scripts/commands/%i_%i.sh'%(slice1,slice2), 'w') as file:
    file.write(filedata)
        
subprocess.run(["sbatch", "/scratch/alpine/mahu8801/blackhole_data/scripts/commands/%i_%i.sh"%(slice1,slice2)],capture_output=True)
    
        