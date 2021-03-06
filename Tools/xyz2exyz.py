import numpy as np

### Read xyz trajectories in <list_names>, convert them in extended xyz for molforge

newf = open("H2O_tr.xyz",'w')

list_names = 'H2O_1_H2O_tr_random H2O_2_H2O_tr_random H2O_3_H2O_tr_random H2O_4_H2O_tr_random H2O_5_H2O_tr_random'

tot_frames = 0
for filename in list_names.split():
    f = open(filename+'.xyz','r')
    lines = f.readlines()
    words = lines[0].split()
    nats = int(words[0])
    nframes = int(len(lines) / (nats+2))
    tot_frames = tot_frames + nframes
    ntypes = 2
    print(filename+'.xyz',nats,nframes)

    box = "100.0 0.0 0.0 0.0 100.0 0.0 0.0 0.0 100.0"

    l = 0
    for i in range(nframes):
        newf.write(lines[l])
        l = l+1
        newf.write(box+" "+str(ntypes)+"\n")
        l = l+1
        for v in range(nats):
            words = lines[l].split()
            if words[0]=="O" :
                newf.write(lines[l].rstrip()+" 1 "+"16.00\n")
                l = l+1
            if words[0]=="H" :
                newf.write(lines[l].rstrip()+" 2 "+"1.000\n")
                l = l+1
    f.close()

newf.close()

print('Total frames:',tot_frames)
