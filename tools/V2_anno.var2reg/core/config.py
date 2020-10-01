'''
Global config

@author: niko
'''

import os, sys, json

if __name__ == '__main__':
    pass

config=json.load(open(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/config/default.tools.json"))

# load user specified config and overloads current tools config.
# Will also replace all configured "prefixes"
def loadUserConfig(confF):
    global config
    user=json.load(open(confF))
    # add tools that were not in user config
    for k in config["tools"].keys():
        if k not in user["tools"].keys():
            user["tools"][k]=config["tools"][k]
    config=user
    if "prefix" in config:
        s=json.dumps(config, separators=(',', ':'))
        for (k, v) in config["prefix"].items():
            s=s.replace(k, v)  
        config=json.loads(s)
    return config

# return pretty print version of current config
def getConfig():
    global config
    return json.dumps(config, indent=4, sort_keys=True)
    
# check for required commandline tools
def checkTools(required):
    global config
    if config is None:
        print("Error: No tools configuration passed! Exiting...")
        sys.exit(1)   

    for t in required :
        if (t not in config["tools"]):
            print("Warn: config for required tool %s was not found! Make sure that %s is a callable commandline..." % t)
    return config

# get tool from config (or return key if not found)
def getTool(key):
    if (key not in config["tools"]):
        return key
    return config["tools"][key]


#loadUserConfig("/Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/config.json")
#print(json.dumps(config, indent=4))