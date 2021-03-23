import os
import time
import sys
import subprocess
from subprocess import Popen, PIPE
# You need `cycles.py` and `to_relax.cfg`

def getCountCfgs():
    command = ' grep "BEGIN_CFG" 4_toRelax/selected.cfg | wc -l '
    if int(os.popen(command).read().split()[0]) >= 1:
        return True
    #
    return False
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

    # command = "cd 5_afterActiveLearning/  && " +\
    #           "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/META/ .  && " +\
    #           "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/saved/ . && " +\
    #           "cd META/  &&  yes | clean.run  &&  cp ../saved/* ."  ## << it doesn't work :/
    command = "cd 5_afterActiveLearning/  && " +\
              "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/META/ .  && " +\
              "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/saved/ . && " +\
              "cd META/  &&  cp ../saved/* ."
    os.system(command)

    command = "cd 4_toRelax/  && " +\
              "cp /home/chinchay/projects/def-rmelnik/chinchay/mydocs/paper1/e_trainedforcarbonwire/4_toRelax/relax.ini .  && " +\
              "cp ../to_relax.cfg ."
    os.system(command)

    command = "cd 2_myTraining/  && " +\
              "cp /home/chinchay/projects/def-rmelnik/chinchay/mydocs/paper1/agnr/1/2_myTraining/pot_blank_binary.mtp .  && " +\
              "cp pot_blank_binary.mtp pot.mtp  && " +\
              "touch train.cfg  && " +\
              "mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg"
    os.system(command)

    command = "cd 4_toRelax/  && " +\
              "cp ../2_myTraining/pot.mtp .  && " +\
              "cp ../2_myTraining/state.mvs .  && " +\
              "mlp relax relax.ini --cfg-filename=to_relax.cfg --min-dist=0.5  && " +\
              "cat selected.cfg_*  > selected.cfg"
    os.system(command)
#

def selectionStep():
    command = "cd 2_myTraining/  && " +\
              "cp ../4_toRelax/selected.cfg .  && " +\
              "mlp select-add pot.mtp train.cfg selected.cfg diff.cfg  && " +\
              "mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar"
    os.system(command)
#

def shouldContinueSleeping():
    command = 'squeue -u chinchay | grep "META" | wc -l'
    nSqueue = int(os.popen(command).read().split()[0])
    if nSqueue >= 1:
        return True
    #
    return False
#

def sendJobs(nJobs):
    command = "cd 5_afterActiveLearning/META/  && " +\
              "submit.run " + str(nJobs)
    os.system(command)
    
    continueSleeping = True
    while continueSleeping:
        time.sleep(120)
        continueSleeping = shouldContinueSleeping()
    #
    ###
    command = "cd 5_afterActiveLearning/META/  && " +\
              "resubmit.run " + str(nJobs)
    string   = "No failed/unfinished jobs; nothing to resubmit"
    st2check = os.popen(command).read().split("\n")[0]
    if st2check != string:
        continueSleeping = True
        while continueSleeping:
            time.sleep(120)
            continueSleeping = shouldContinueSleeping()
        #    
    ######

    # The following did not work because I cannot get the correct signal for waiting until finished jobs
    # I tried to implement the advice here https://blog.miguelgrinberg.com/post/how-to-make-python-wait
    # other advices here https://hpc-discourse.usc.edu/t/how-to-wait-until-the-job-is-finished/242/2
    # command = 'cd 5_afterActiveLearning/META/  && ' +\
    #           'echo "submit.run ' + str(2) + ' --wait &" > myexecute.sh  && ' +\
    #           'echo "wait" >> myexecute.sh  && ' +\
    #           'bash myexecute.sh'
    # os.system(command)

    # ########
    # # wait

    # command = 'cd 5_afterActiveLearning/META/  && ' +\
    #           'echo "resubmit.run ' + str(2) + ' --wait &" > myexecute.sh  && ' +\
    #           'echo "wait" >> myexecute.sh  && ' +\
    #           'bash myexecute.sh'
    # os.system(command)

    ########
    # wait
    # if "nothing to resubmit" continue
#

def dftStep(nJobs):
    command = "cd 5_afterActiveLearning  && rm -rf META/"
    os.system(command)

    command = "cd 5_afterActiveLearning/  && " +\
              "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/META/ .  && " +\
              "cd META/  &&  cp ../saved/* .  && " +\
              "cp ../../2_myTraining/POSCAR* . "
    os.system(command)

    command = "ls 5_afterActiveLearning/META/POSCAR* | wc -l"
    nPoscars = int(os.popen(command).read().split()[0])

    os.system("rm -f 5_afterActiveLearning/META/table.dat")
    f = open("5_afterActiveLearning/META/table.dat", "a")
    for i in range(nPoscars):
        command = "cd 5_afterActiveLearning/META/  &&  python vasp2qe.py -i \"POSCAR" + str(i) + "\""
        os.system(command)
        
        line = "srun pw.x < ../scf_" + str(i + 1) + ".in > scf_" + str(i + 1) + ".out\n"
        f.write(line)
    #
    f.close()
    #
    print("sending jobs now...")
    sendJobs(nJobs)
    print("It seems I have finished jobs")
#

def wait2SeeIfTrainGotStuck(lastLine0, lastLine, checkTrainTime=120):
    endTrainingLine = "_______________________________________________"   
    command = "tail -2 training.txt"
    # First, check if trainig has finished, so it is already unstuck.
    if endTrainingLine == os.popen(command).read().split()[0]:
        return False 
    #
    # sleep and then compare lines:
    time.sleep(checkTrainTime)
    return (lastLine0 != lastLine)
#

def isTrainingUnstuck(myprocess, checkTrainTime):
    # https://stackoverflow.com/questions/12057794/python-using-popen-poll-on-background-process
    isUnstuck = True
    lastLine0 = ""
    # `myprocess.poll()` will check if subprocess is finished. If not, it will return `None`
    while (myprocess.poll() is None) and isUnstuck: # BOTH options MUST be true
        lastLine  = os.popen("tail -1 training.txt").read().split()[0]
        isUnstuck = wait2SeeIfTrainGotStuck(lastLine0, lastLine, checkTrainTime)
        lastLine0 = lastLine
    #
    return isUnstuck 
#

def hasTrainingFinishedAndUnstuck(checkTrainTime):
    f = open("training.txt", "w")
    g = open("errorsByPythonPopen.txt", "w")
    p = Popen(["cd", "2_myTraining/", "&&", "mlp", "train", "pot.mtp", "train.cfg"], stdout=f, stderr=g)
    #
    unStuck = isTrainingUnstuck(p, checkTrainTime) # <<== don't worry, it sleeps until finding Trainig has finished, or it got stuck
    if not unStuck:  # it got stuck, so It couldn't finish calculations
        p.kill()
    #
    f.close()
    g.close()
    return unStuck
#

def checkErrorFile():
    g = open("errorsByPythonPopen.txt", "r")
    checkLine = g.readline().split(" \n")[0]
    g.close()
    return checkLine
#

def trainingStep(checkTrainTime):
    command = "cd 5_afterActiveLearning/META/  && " +\
              "python scriptQE2cfg.py  && " +\
              "cd ../../2_myTraining/  && " +\
              "rm -f POSCAR* && rm -f diff.cfg && rm -f temp1.cfg && rm -f selected.cfg && rm -f training.txt  && " +\
              "cp ../5_afterActiveLearning/META/train.cfg train2.cfg  && " +\
              "cat train2.cfg >> train.cfg  && " +\
              "rm train2.cfg"
    os.system(command)
    #
    isTrainingOK = False
    for _ in range(2):        
        if not isTrainingOK:
            if hasTrainingFinishedAndUnstuck(checkTrainTime):
                # mlp exited, and produced a `training.txt`. However we need to check if there was a BGFS ascending error:
                isTrainingOK = ( checkErrorFile() != errorBFGSaccending )
            else:
                # mlp got stuck, runs infinitely. `training.txt` has stuck in the same line
                isTrainingOK = False
            #
        #
    #
    if not isTrainingOK:
        # let's try another trick:
        # if the same error appear again, take a fresh pot.mtp:
        command = "cp pot_blank_binary.mtp pot.mtp"
        os.system(command)
        if hasTrainingFinishedAndUnstuck(checkTrainTime):
            if checkErrorFile() == errorBFGSaccending:
                #? You should increase configs in 4_toRelax/to_reala.cfg, and start all again
                print("error errorBFGSaccending again?")
                print("You should increase configs in 4_toRelax/to_reala.cfg, and start all again")
                print("I am going to stop here... :/")
                sys.exit()
            #
        #
    #
    #
    command = "cd 2_myTraining/  &&  mv Trained.mtp_ pot.mtp  && " +\
              "mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg"
    os.system(command)
#

def relaxStep():
    command = "cd 4_toRelax/  && " +\
              "rm -f select*  && " +\
              "cp ../2_myTraining/pot.mtp .  && " +\
              "cp ../2_myTraining/state.mvs .  && " +\
              "mlp relax relax.ini --cfg-filename=to_relax.cfg --min-dist=0.5  && " +\
              "cat selected.cfg_*  > selected.cfg"
    os.system(command)
#


# os.system("mlj")##??

#initialize()
#
maxNcycles = 1
nJobs = 1
checkTrainTime = 120
continuar = True

for i in range(maxNcycles):
    if continuar:
        #selectionStep()
        #dftStep(nJobs)
        trainingStep(checkTrainTime)
        relaxStep()
        continuar = shouldContinue()
        print("end of loop i = " + str(i) )
#     #
print("finish loop")
# #

