#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 16:11:19 2019

@author: dobri
"""

import xwiimote
from select import poll, POLLIN
import time

def display_some_info_and_open_iface(dev,num):
    try:
        fd = dev.get_fd()
        dev.open(dev.available() | xwiimote.IFACE_WRITABLE)
        mode_num_temp = num
        for ledn in range(0,4):
            dev.set_led(ledn+1,bool(mode_num_temp%2))
            mode_num_temp=(mode_num_temp//2)
    
        print "Capacity", '\t\t',  dev.get_battery(), "%"
        print '\n'
        
    except SystemError as e:
        print "ooops"
        e
        exit(1)
    return fd

def calibrate_acc(dev,pwiimote):
    global OFFSET_X
    global OFFSET_Y
    global OFFSET_Z
    global OFFSET_X_3
    global OFFSET_Y_2
    global OFFSET_Z_1
    print "syspath:" + dev.get_syspath()
    print "opened mask:", dev.opened()
    print "opened mask:", dev.opened()
    print "devtype:",   dev.get_devtype()
    print "extension:", dev.get_extension()
    print "Rest the wiimote face-up"
    dev.set_led(2, True)
    rest_time = time.time()    
    while time.time() - rest_time < 10:
        polls = pwiimote.poll()
        for fd, evt in polls:
            if fd == dev.get_fd():
                dev.dispatch(revt)
                if revt.type == xwiimote.EVENT_ACCEL:
                    x1,y1,z1 = revt.get_abs(0)

    dev.set_led(2, False)    
    print "Rest the wiimote butt-up"
    dev.set_led(3, True)
    rest_time = time.time()    
    while time.time() - rest_time < 10:
        polls = pwiimote.poll()
        for fd, evt in polls:
            if fd == dev.get_fd():
                dev.dispatch(revt)
                if revt.type == xwiimote.EVENT_ACCEL:
                    x2,y2,z2 = revt.get_abs(0)

    dev.set_led(3, False)    
    print "Rest the wiimote on the wiimote's right shoulder"
    dev.set_led(4, True)
    rest_time = time.time()    
    while time.time() - rest_time < 10:
        polls = pwiimote.poll()
        for fd, evt in polls:
            if fd == dev.get_fd():
                dev.dispatch(revt)
                if revt.type == xwiimote.EVENT_ACCEL:
                    x3,y3,z3 = revt.get_abs(0)
    
    OFFSET_X_3 = x3
    OFFSET_Y_2 = y2
    OFFSET_Z_1 = z1
    OFFSET_X = float((x1 + x2) / 2)
    OFFSET_Y = float((y1 + y3) / 2)
    OFFSET_Z = float((z2 + z3) / 2)
    dev.set_led(4, False)    
    print "Done!"
    print "The new offsets in x, y, and z are ", OFFSET_X, OFFSET_Y, OFFSET_Z, "."
    print "The other supporting values of x, y, and z are ", OFFSET_X_3, OFFSET_Y_2, OFFSET_Z_1, "."
##
    
    
    
if __name__ == '__main__':
    
    # For some reason it makes an empty slot, which then confuses everything.
    # So hardcode this.
    if True:
        subtract_one_weird_wii_index = 1
    else:
        subtract_one_weird_wii_index = 0
    
    
    prog_start = time.time()

    print "=== " + xwiimote.NAME_CORE + " ==="
    
    # List and remember wiimotes
    try:
        monitor = xwiimote.monitor(True, True)
        wiis = list()
        ent = monitor.poll()
        wiis.append(ent)
        while ent is not None:
            print "Found device: " + ent
            ent = monitor.poll()
            wiis.append(ent)
    except SystemError as e:
            print "oooooops cannot create monitor (", e, ")"
    
    
    # Continue only if there are two wiimotes
    #print len(wiis)
    #print wiis[0]
    #print wiis[1]
    #wiis=wiis[0]
    #print wiis
    #print len(wiis)
    # Which of the found to use?
    #d=0
    if len(wiis) < 2:
        print "Less than two wiimotes to read!"
        exit(0)
    
    # Create two new ifaces
    try:
        devs = list()
        for d in range(0,len(wiis)-subtract_one_weird_wii_index):
            dev = xwiimote.iface(wiis[d])
            devs.append(dev)
    except IOError as e:
        print "ooops,"
        e
        exit(1)
    
    
    # Display some information and open the iface
    fds = list()
    for d in range(0,len(wiis)-subtract_one_weird_wii_index):
        f = display_some_info_and_open_iface(devs[d],d+1)
        fds.append(f)
    
        
    # read some values
    pwiimote = poll()
    for d in range(0,len(wiis)-subtract_one_weird_wii_index):
        pwiimote.register(devs[d].get_fd(), POLLIN)
    
    revt = xwiimote.event()
    
    for d in range(0,len(wiis)-subtract_one_weird_wii_index):
        print wiis[d]
        calibrate_acc(devs[d],pwiimote)
    
    if 0:
        #rest_time = time.time()    
        #while time.time() - rest_time < 5:
        try:
            while True:
                polls = pwiimote.poll()
                for fd, evt in polls:
                    if fd == fds[0]:
                        devs[0].dispatch(revt)
                        if revt.type == xwiimote.EVENT_ACCEL:
                            x1,y1,z1 = revt.get_abs(0)
                            print fd, "ACC", "%8.2f" % x1, "%8.2f" % y1, "%8.2f" % z1
                    if fd == fds[1]:
                        devs[1].dispatch(revt)
                        if revt.type == xwiimote.EVENT_ACCEL:
                            x2,y2,z2 = revt.get_abs(0)
                            print '\t\t\t\t\t\t', fd, "ACC", "%8.2f" % x2, "%8.2f" % y2, "%8.2f" % z2
        except KeyboardInterrupt:
            print "exiting..."

    # cleaning
    for d in range(0,len(fds)):
        pwiimote.unregister(fds[d])
    #exit(0)