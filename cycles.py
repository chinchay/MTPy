import os
import time
import sys
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

def trainingStep():
    command = "cd 5_afterActiveLearning/META/  && " +\
              "python scriptQE2cfg.py  && " +\
              "cd ../../2_myTraining/  && " +\
              "rm -f POSCAR* && rm -f diff.cfg && rm -f temp1.cfg && rm -f selected.cfg && rm -f training.txt  && " +\
              "cp ../5_afterActiveLearning/META/train.cfg train2.cfg  && " +\
              "cat train2.cfg >> train.cfg  && " +\
              "rm train2.cfg"
    os.system(command)
    
    errorBFGSaccending = "ERROR: BFGS: stepping in accend direction detected."
    command = "cd 2_myTraining/  &&  mlp train pot.mtp train.cfg > training.txt"
    if errorBFGSaccending == os.popen(command).read().split("\n")[0]:
        # try again:
        command = "mlp train pot.mtp train.cfg > training.txt"
        # if the same error appear again, take a fresh pot.mtp:
        if errorBFGSaccending == os.popen(command).read().split("\n")[0]:
            command = "cp pot_blank_binary.mtp pot.mtp  && " +\
                      "mlp train pot.mtp train.cfg > training.txt"
            #? error again? You should increase configs in 4_toRelax/to_reala.cfg, and start all again
            if errorBFGSaccending == os.popen(command).read().split("\n")[0]:
                print("error errorBFGSaccending again?")
                print("You should increase configs in 4_toRelax/to_reala.cfg, and start all again")
                print("I am going to stop here... :/")
                sys.exit()
            #
        #
    # elif
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
continuar = True

for i in range(maxNcycles):
    if continuar:
        #selectionStep()
        #dftStep(nJobs)
        trainingStep()
        relaxStep()
        continuar = shouldContinue()
        print("end of loop i = " + str(i) )
#     #
print("finish loop")
# #

