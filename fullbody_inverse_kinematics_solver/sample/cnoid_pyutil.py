from cnoid.Base import *
from cnoid.BodyPlugin import *
#from cnoid.PythonSimScriptPlugin import *
import cnoid.Body

import numpy as np

def getItemTreeView():
    if callable(ItemTreeView.instance):
        return ItemTreeView.instance()
    else:
        return ItemTreeView.instance


def getRootItem():
    if callable(RootItem.instance):
        return RootItem.instance()
    else:
        return RootItem.instance


def getWorld():
    rI = getRootItem()
    ret = rI.findItem("World")
    if ret == None:
        ret = WorldItem()
        rI.addChildItem(ret)
        getItemTreeView().checkItem(ret)
    return ret


def cnoidPosition(rotation = None, translation = None):
  ret = np.identity(4)
  if not (rotation is None):
    ret[:3, :3] = rotation
  if not (translation is None):
    ret[:3, 3] = translation
  return ret


def cnoidRotation(cPosition):
  return cPosition[:3, :3]


def cnoidTranslation(cPosition):
  return cPosition[:3, 3]


def loadRobot(fname, name = None):
    print('loadRobot: %s'%(fname))
    bI = BodyItem()
    bI.load(fname)
    if name:
        bI.setName(name)
    if callable(bI.body):
        rr = bI.body()
    else:
        rr = bI.body
    rr.calcForwardKinematics()
    bI.storeInitialState()
    wd = getWorld()
    if callable(wd.childItem):
        wd.insertChildItem(bI, wd.childItem())
    else:
        wd.insertChildItem(bI, wd.childItem)
    getItemTreeView().checkItem(bI)
    return rr


#def loadRobot(fname, name = None):
#    print('loadRobot: %s'%(fname))
#    import time
#    time.sleep(0.5)
#    return fname

def findItem(name):
    return getRootItem().findItem(name)

def findRobot(name):
    ret = findItem(name)
    ## add class check...
    if callable(ret.body):
        return ret.body()
    else:
        return ret.body

def flushRobotView(name):
    findItem(name).notifyKinematicStateChange()
    #MessageView.getInstance().flush()
    MessageView.instance.flush()

