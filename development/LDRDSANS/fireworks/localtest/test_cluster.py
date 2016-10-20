from fireworks import Firework, LaunchPad, ScriptTask
from fireworks.core.rocket_launcher import launch_rocket, rapidfire

# set up the LaunchPad and reset it
launchpad = LaunchPad(strm_lvl='WARNING') # set messaging lowest level to WARNING
launchpad.reset('', require_password=False)

# create the Firework consisting of a single task
firetask = ScriptTask.from_str('cd /projects/development/LDRDSANS/fireworks/localtest; ./test_cluster.sh')
firework = Firework(firetask)
fw_yaml = firework.to_file("my_firework.yaml") # save to yaml file, and get the string representation
fw_json = firework.to_file("my_firework.json") # save to json file, and get the string representation

# store workflow and launch it locally
launchpad.add_wf(firework)
launch_rocket(launchpad)  # same as "rlaunch singleshot"
#rapidfire(launchpad, FWorker(), strm_lvl='WARNING')  # same as "rlaunch rapidfire"

# loading from file
# any class in FireWorks that subclasses FWSerializable has methods from_file and to_file
#firework = Firework.from_file("fw_test.yaml")
#fw_yaml = Firework.from_file("fw_test.json")
