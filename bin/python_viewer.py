#!/usr/bin/python

from visual import *

scene.autoscale=False
scene.autocenter=False
sphere()
down=False
while 1:
    rate(100)
    if scene.mouse.events:
        m=scene.mouse.getevent()
        if m.press=='left':
            down=True
            lastpos=scene.mouse.pos
        elif m.release=='left':
            down=False
    if down and scene.mouse.pos!=lastpos:
        print(scene.mouse.pos)
        scene.center=scene.center+(scene.mouse.pos-lastpos)
        lastpos=2*scene.mouse.pos-lastpos
        
