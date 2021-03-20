import os

def getCountCfgs():
    command = "count selected.cfg"
    return os.popen(command).readlines().split()[0]
#

def shouldContinue():
    if getCountCfgs() == 0:
        return False
    #
    return True
#

def initialize():
    os.system("mkdir 2_myTraining/")
    os.system("mkdir 4_toRelax/")
    os.system("mkdir 5_afterActiveLearning/")

    os.system("cd 5_afterActiveLearning/")
    os.system("cp -r ../../../META/ .") ###################<<<<<
    os.system("cp -r ../../b/5_afterActiveLearning/saved/ .") ##<<<<<<<<
    os.system("cd META/")
    os.system("yes | clean.run")
    os.system("cp ../saved/* .")

    os.system("cd ../../2_myTraining/")
    os.system("cp ../../b/4_toRelax/beforeRelax1/pot_blank.mtp pot.mtp") #<<<<<<<
    os.system("cp ../../b/4_toRelax/relax.ini 4_toRelax/relax.in") #<<<<<<


    os.system("touch train.cfg")
    os.system("mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg")

    os.system("cd ../4_toRelax/")
    os.system("cp ../2_myTraining/pot.mtp .")
    os.system("cp ../2_myTraining/state.mvs .")
    os.system("mlp relax relax.ini --cfg-filename=to_relax.cfg --min-dist=0.5")
    os.system("cat selected.cfg_*  > selected.cfg")
#
def selectionStep():
    os.system("cd ../2_myTraining/")
    os.system("cp ../4_toRelax/selected.cfg .")
    os.system("mlp select-add pot.mtp train.cfg selected.cfg diff.cfg")
    os.system("mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar")
#

def dftStep():
    os.system("cd ../5_afterActiveLearning/META/")
    os.system("yes | clean.run")
    os.system("cp ../saved/* .")
    os.system("cp ../../2_myTraining/POSCAR* .")

    os.system("rm table.dat")
    f = open("table.dat", "a")
    for i in range(n):
        command = "python vasp2qe.py -i \"POSCAR" + str(i) + "\""
        os.system(command)
        
        line = "srun pw.x < ../scf_" + str(i) + ".in > scf_" + str(i) + ".out\n"
        f.write(line)
    #
    f.close()
    #
    sendJobs()
#

def trainingStep():
    os.system("python scriptQE2cfg.py")
    os.system("cd ../../2_myTraining/")
    os.system("rm POSCAR* && rm diff.cfg && rm temp1.cfg")
    os.system("cp ../5_afterActiveLearning/META/train.cfg train2.cfg")
    os.system("cat train2.cfg >> train.cfg")
    os.system("rm train2.cfg")
    os.system("mlp train pot.mtp train.cfg > training.txt")
    os.system("mv Trained.mtp_ pot.mtp")
    os.system("mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg")
#

def relaxStep():
    os.system("cd ../4_toRelax/")
    os.system("cp ../2_myTraining/pot.mtp .")
    os.system("cp ../2_myTraining/state.mvs .")
    os.system("mlp relax relax.ini --cfg-filename=to_relax.cfg --min-dist=0.5")
    os.system("cat selected.cfg_*  > selected.cfg")
#

def sendJobs():
    command ="submit.run " + str(10) + " &  && wait" #<<<<<<
    os.system(command)

    ########
    # wait

    command = "resubmit.run " + str(10) + " &  && wait" #<<<<<<
    os.system(command)

    ########
    # wait
    # if "nothing to resubmit" continue
#



initialize()
#
maxNcycles = 2
continuar = True
for _ in range(maxNcycles):
    if continuar:
        selectionStep()
        dftStep()
        trainingStep()
        relaxStep()
        continuar = shouldContinue()
    #
#

