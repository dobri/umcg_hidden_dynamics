# -*- coding: utf-8 -*-
"""
@Author: Dobri Dotov, @C3, UNAM, 2016
"""
import numpy as np
import pyaudio
import math
import time
import wave
from select import poll, POLLIN
import sys, getopt
from pymouse import PyMouse
import os


def map_eff_states_to_screen_bars(x0,x1):
    BARTHETA[0] = x0
    BARTHETA[1] = x1 #(x1-(PI/2))*1+(PI/2)
    if BARTHETA[0]>PI:
        BARTHETA[0] = (BARTHETA[0] % PI)*-1 + PI
    BARTHETA[0] = (BARTHETA[0]-PI/2)*SCALE_KUR_AMP+(PI/2)
    #if BARTHETA[1]>PI:
    #    BARTHETA[1] = (BARTHETA[1] % PI)*-1 + PI
    return BARTHETA


def convert_to_angle(s,k):
    # The participant input as angle.
    theta = s[k-1]+PI
    if theta>PI:
        theta=PI
    if theta<0:
        theta=0
    return theta

        
def convert_to_angle_2(s,k):
    if k>19:
        if arduino_status:
            theta_smoothed[0] = np.mean(s[(k-19):(k-1)])
            theta_smoothed[1] = np.mean(s[(k-20):(k-2)])
        else:
            theta_smoothed[0] = np.mean(s[(k-3):(k-1)])
            theta_smoothed[1] = np.mean(s[(k-4):(k-2)])
    theta = np.arctan2(np.diff(theta_smoothed)/dt_ceil,s[k-1])
    if theta<0:
        theta = (2*PI)-theta
    return theta


def compute_and_return_performance_feedback(log_file_name,col1name,col2name):
    from wcc_by_rmse import wcc_by_rmse
    import pandas as pd
    X = pd.read_csv(log_file_name,skipinitialspace=True)
    fps=1/np.mean(np.diff(X['Time']))
    score = 1e2*wcc_by_rmse(X[col1name][int(fps*5):(len(X)-1)],X[col2name][int(fps*5):(len(X)-1)],X['Time'][int(fps*5):(len(X)-1)],1,5,False)
    try:
        contents=pd.read_csv('scores',header=None,na_values='   nan')
    except:
        contents=[]
    f = open('scores','a+')
    f.write("%6.3f\n" % score) # "%s," % log_file_name, 
    f.close()
    if len(contents)>0:
        print "%s\n" % 'Your scores from previous trials:'
        for c in contents.values:
            try:
                print "%6.3f" % c
            except:
                print "%6s" % c
    print "%s\n" % 'Your score from the last trial:'
    print "%6.3f" % score

    import matplotlib.pyplot as plt
    font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 20}
    plt.rc('font', **font)
    plt.rc('axes', labelsize=30)

    fig  = plt.figure(figsize=(10,10))
    axes = fig.add_axes([.15,.15,.75,.75])
    if len(contents)>0:
        axes.plot(contents,'-o',linewidth=2,markersize=15,label='Previous Trials')
    axes.plot(len(contents),np.array(score),'-o',linewidth=2,markersize=15,color='r',label='Last Trial')
    axes.set_xlabel('Trial')
    axes.set_xticks(range(0,len(contents)+1,1))
    axes.set_xticklabels(range(1,len(contents)+2,1))
    axes.set_ylabel('Score')
    #axes.set_yticks(range(np.floor(contents.min()),np.ceil(contents.max())+1,10))
    #axes.set_ylim(np.floor(contents.min())-10,np.ceil(contents.max())+10)
    axes.set_title('Performance (Synchronization and Pitch-matching)')
    axes.legend(loc='upper center', shadow=False, fontsize='medium')
    plt.show()
    
    
def print_sound_device_info(p):
    max_apis = p.get_host_api_count()
    max_devs = p.get_device_count()
  
    print("PortAudio System Info:")
    print("=" * 120)
    print("Version: %d" % pyaudio.get_portaudio_version())
    print("Version Text: %s" % pyaudio.get_portaudio_version_text())
    print("Number of Host APIs: %d" % max_apis)
    print("Number of Devices  : %d" % max_devs)
  
    print("\nHost APIs:\n==========")
  
    for i in range(max_apis):
        apiinfo = p.get_host_api_info_by_index(i)
        for k in apiinfo.items():
            print("%s: %s" % k)
        print("--------------------------")
  
    devinfo = p.get_device_info_by_index(0)
    NUM_CHANS=devinfo['maxOutputChannels']
    print(NUM_CHANS)
    
    print("\nDevices:\n========")
  
    for i in range(max_devs):
        devinfo = p.get_device_info_by_index(i)
        print(devinfo['maxOutputChannels'])
        
        # print out device parameters
        for k in devinfo.items():
            name, value = k
            # if host API, then get friendly name
            if name == 'hostApi':
                value = str(value) + \
                    " (%s)" % p.get_host_api_info_by_index(k[1])['name']
            print("\t%s: %s" % (name, value))
    return NUM_CHANS

def sigmoid_fun(sigm_in,sigm_a,sigm_b,sigm_c,sigm_d,sigm_e):
    return sigm_a + sigm_b/(1 + math.exp( - sigm_c * ( sigm_d * sigm_in - sigm_e )))

def rescale_x_acc_fun_linear(x_acc_in):
    if math.fabs(x_acc_in) <= 1:
        x_acc_in = math.asin(x_acc_in/G)/(PI/2)
    else:
        x_acc_in = x_acc_in 
    return x_acc_in

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
        print "ooops", e
        exit(1)
    return fd

def rescale_acc(x,d=0):
    return (x - OFFSETS[0,:,d])/(OFFSETS[1,:,d]-OFFSETS[0,:,d])

def rescale_ard_acc(x):
    return (x - OFFSETS_ARD[0,:])/(OFFSETS_ARD[1,:]-OFFSETS_ARD[0,:])
    
def set_XL_to(setx):
    global XL
    XL = setx
    
def set_VL_to(setv):
    global VL
    VL = setv

def set_XR_to(setx):
    global XR
    XR = setx
    
def set_VR_to(setv): #,setv2,setv3
    global VR #,VR_2,VR_3
    VR = setv

def set_XL2_to(setx):
    global XL2
    XL2 = setx

def set_VL2_to(setx):
    global VL2
    VL2 = setx

def set_XR2_to(setx):
    global XR2
    XR2 = setx

def set_VR2_to(setx):
    global VR2
    VR2 = setx

def set_BUFFERS_N(bn):
    global BUFFERS_N
    BUFFERS_N = bn

def set_FREQL0(fr):
    global FREQL0
    FREQL0 = fr
            
def set_FREQR0(fr):
    global FREQR0
    FREQR0 = fr

def set_FREQL02(fr):
    global FREQL02
    FREQL02 = fr

def set_FREQR02(fr):
    global FREQR02
    FREQR02 = fr

def set_XX(xvec):
    global XX
    XX[((BUFFERS_N - 1) * CHUNK + 0):((BUFFERS_N) * CHUNK - 1)] = xvec

def cdo(xin,yin,DT,freq):
    xout=xin+DT*yin
    yout=yin+DT*(-(freq*2*math.pi)**2*xin-GAMMA*yin*((freq*2*math.pi)**2*xin**2/2+yin**2/2-BETA))
    return (xout,yout)

def osc():
    global WAVEDATAL
    global WAVEDATAR
    global XL,VL,XR,VR
    for n in range(0, CHUNK, 1):
        FREQVECL[n]=FREQL0 + (n / (CHUNK - 1.00)) * (FREQL1 - FREQL0)
        FREQVECR[n]=FREQR0 + (n / (CHUNK - 1.00)) * (FREQR1 - FREQR0)
        
        if AMPLIFY_OSC_L>0:
            xl=XL
            vl=VL
            for h in range(0,len(Fs_RATIOS[0]),1):
                if XL[h]>1:
                    XL[h]=.1
                if VL[h]>100:
                    VL[h]=.1
                xl[h]=XL[h] + DT*VL[h]
                vl[h]=VL[h] + DT*(-(FREQVECL[n]*Fs_RATIOS[0,h]*2*math.pi)**2*XL[h]-GAMMA*VL[h]*((FREQVECL[n]*Fs_RATIOS[0,h]*2*math.pi)**2*XL[h]**2/2+VL[h]**2/2-BETA))
            WAVEDATAL[n] = sum(xl)/(h+1)
            XL=xl
            VL=vl

        if AMPLIFY_OSC_R>0:
            xr=XR
            vr=VR
            for h in range(0,len(Fs_RATIOS[1]),1):
                fr = min((3000,FREQVECR[n]*Fs_RATIOS[1,h]))
                if XR[h]>1:
                    XR[h]=.1
                if VR[h]>100:
                    VR[h]=.1
                xr[h]=XR[h] + DT*VR[h]
                vr[h]=VR[h] + DT*(-(fr*2*math.pi)**2*XR[h]-GAMMA*VR[h]*((fr*2*math.pi)**2*XR[h]**2/2+VR[h]**2/2-BETA))
            WAVEDATAR[n] = sum(xr)/(h+1)
            XR=xr
            VR=vr
        
    set_FREQL0(FREQL1)
    set_FREQR0(FREQR1)

    if NUM_CHANS>2:
        global WAVEDATAL2
        global WAVEDATAR2
        XL02=XL2
        VL02=VL2
        XR02=XR2
        VR02=VR2
        for n in range(0, CHUNK, 1):
            if AMPLIFY_OSC_L>0:
                FREQVECL2[n]=FREQL02 + (n / (CHUNK - 1.00)) * (FREQL12 - FREQL02)
                xl=XL02 + DT*VL02
                vl=VL02 + DT * (-(FREQVECL2[n] * 2 * math.pi) ** 2 * XL02 - GAMMA * VL02 * ((FREQVECL2[n] * 2 * math.pi) ** 2 * XL02 ** 2 / 2 + VL02 ** 2 / 2 - BETA))
                WAVEDATAL2[n] = xl
                XL02=xl
                VL02=vl

            if AMPLIFY_OSC_R>0:
                FREQVECR2[n]=FREQR02 + (n / (CHUNK - 1.00)) * (FREQR12 - FREQR02)
                xr=XR02 + DT*VR02
                vr=VR02 + DT * (-(FREQVECR2[n] * 2 * math.pi) ** 2 * XR02 - GAMMA * VR02 * ((FREQVECR2[n] * 2 * math.pi) ** 2 * XR02 ** 2 / 2 + VR02 ** 2 / 2 - BETA))
                WAVEDATAR2[n] = xr
                XR02=xr
                VR02=vr

        set_FREQL02(FREQL12)
        set_XL2_to(xl)
        set_VL2_to(vl)
        set_FREQR02(FREQR12)
        set_XR2_to(xr)
        set_VR2_to(vr)


def callback2(in_data, frame_count, time_info, status):
    set_BUFFERS_N(BUFFERS_N + 1)
    osc() # This sets the WAVEDATA sound buffer using global params.
    chunk  = np.ravel(np.vstack((np.multiply(WAVEDATAL,AMPLIFY_OSC_L), np.multiply(WAVEDATAR,AMPLIFY_OSC_R))), order='F')
    if SOUND_EXPORT_FLAG:
        global sound_frames
        sound_frames.append(chunk)
    chunk = chunk.astype(np.float32).tostring()
    return (chunk, pyaudio.paContinue)


def callback4(in_data, frame_count, time_info, status):
    set_BUFFERS_N(BUFFERS_N + 1)
    osc() # This sets the WAVEDATA sound buffer using global params.
    chunk  = np.ravel(np.vstack((np.multiply(WAVEDATAL,AMPLIFY_OSC_L), np.multiply(WAVEDATAR,AMPLIFY_OSC_R), np.multiply(WAVEDATAL2,AMPLIFY_OSC2), np.multiply(WAVEDATAR2,AMPLIFY_OSC2))), order='F')
    if SOUND_EXPORT_FLAG:
        global sound_frames
        sound_frames.append(chunk)
    chunk = chunk.astype(np.float32).tostring()
    return (chunk, pyaudio.paContinue)


def CartPole(dt,XvecCart,XvecCOM,force_gain=1.,limitlow=-.9,limitup=.9):
    #print(XvecCOM,XvecCart)
    XvecCOM = XvecCOM + dt*(XvecCOM-XvecCart)*force_gain
    # Contain the x-axis range
    XvecCOM = max(limitlow,XvecCOM)
    XvecCOM = min(limitup,XvecCOM)
    return XvecCOM


def Lorenz(dt,XvecLorenz,ForcingLorenz,rho=28,sigma=10,beta=8/3):
    XvecLorenz[0]=XvecLorenz[0]+dt*(sigma*(XvecLorenz[1]-XvecLorenz[0]))
    XvecLorenz[1]=XvecLorenz[1]+dt*(  rho*XvecLorenz[0]-XvecLorenz[1]-XvecLorenz[0]*XvecLorenz[2])
    XvecLorenz[2]=XvecLorenz[2]+dt*(-beta*XvecLorenz[2]+XvecLorenz[0]*XvecLorenz[1]+ForcingLorenz)
    return XvecLorenz


def Chua1(dt,XvecChua,ForcingChua,g=-1,g0=1,b1=9,b2 =5,b3=1):
    g0=g0 - ForcingChua
    #print XvecChua
    XvecChua[0]=XvecChua[0]+dt*( -g *XvecChua[0] - b1*( XvecChua[0] + XvecChua[2] ) * XvecChua[1]**2 )
    XvecChua[1]=XvecChua[1]+dt*( -g0*XvecChua[1] + b2*( XvecChua[0]**2 - XvecChua[2]**2) * XvecChua[1] )
    XvecChua[2]=XvecChua[2]+dt*( -g0*XvecChua[2] + b3*( XvecChua[0] + XvecChua[2] ) * XvecChua[1]**2 )
    
    # Here is where the break happens?
    if XvecChua[1]<0:
        XvecChua[1]=XvecChua[1]*-1
    return XvecChua   


def Chua2(dt,XvecChua,ForcingChua,alphaC=15.6,betaC=28,m0=-1.143,m1=-.714):
    h = m1 * XvecChua[0] + .5*(m0 - m1)*(math.fabs(XvecChua[0] + 1) - math.fabs(XvecChua[0]-1))
    XvecChua[0]=XvecChua[0]+dt*( alphaC*(XvecChua[1] - XvecChua[0] - h) )
    XvecChua[1]=XvecChua[1]+dt*( XvecChua[0] - XvecChua[1] + XvecChua[2] )
    XvecChua[2]=XvecChua[2]+dt*( -betaC*XvecChua[1] )
    return XvecChua   


def artificial_SensorGain_sine(artif_omega,artif_t):
    return (math.cos(artif_omega*artif_t)+1)/3


def midi_key_to_hz(key):
    output_hz = 2**((key - 69) / 12)*440
    
    # Is it here where the singularity happens?
    if output_hz>6000:
        print "Attention! Sound clipping."
        output_hz>6000
    if output_hz<50:
        print "Attention! Sound clipping."
        output_hz=50
    return output_hz


def map_x_to_pixel(x):
    return np.multiply(x,vis_mod_center[0])


def map_x_to_note(xtonote):
    return xtonote*5 + 60


def first_camera_frame(cap):
    global hsv
    global xrange0, xrange1
    global of_fb_winsize
    global center
    
    cap.set(4,240*1) # 1080 720
    cap.set(3,320*1) # 1920 1280

    ret,frame0 = cap.read()
    s = np.shape(frame0)
    hsv = np.zeros_like(frame0)
    hsv[...,1] = 255
    time.sleep(.1)
    
    ret,frame0 = cap.read()
    if GPU_FLAG:
        global cuframe0
        cuframe0 = cv2.cuda_GpuMat()
        cuframe1 = cv2.cuda_GpuMat()
        cuframe0.upload(frame0)
        cuframe0 = cv2.cuda.cvtColor(cuframe0,cv2.COLOR_BGR2GRAY)
        cuframe0 = cv2.cuda.flip(cuframe0,1) # It turns out we have to flip the array.
    
    frame0 = cv2.cvtColor(frame0,cv2.COLOR_BGR2GRAY)
    frame0 = cv2.flip(frame0,1) # It turns out we have to flip the array.
    
    print "Your video frame size is %d by %d." % s[0:2]
    of_fb_winsize = np.mean(np.divide(s[0:2],30),dtype='int')
    center=(np.int(np.round(s[1]/2)),np.int(np.round(s[0]/2)))
    
    
    if GPU_FLAG:
        # Prepare the optic flow computer
        global cuda_farn
        cuda_farn = cv2.cuda.FarnebackOpticalFlow_create(5,.5,False,of_fb_winsize,3,5,1.1,0)
    
    
    # Define a hard-coded ROI. Here it's just left and right sides of the frame.
    #xrange0 = list(range(0,np.int(np.round(s[1]/2))))
    #xrange1 = list(range(np.int(np.round(s[1]/2)),s[1]-1))
    xrange0 = list(range(0,s[1]-1))
    xrange1 = list(range(s[1]-1,s[1]-1))
    
    # Allow the world to move a little bit relative to the previous snapshot
    # and get a new snapshot.
    time.sleep(.1)
    ret,frame1 = cap.read()
    
    if GPU_FLAG:
        global flow_gpu
        cuframe1.upload(frame1)
        cuframe1 = cv2.cuda.cvtColor(cuframe1,cv2.COLOR_BGR2GRAY)
        cuframe1 = cv2.cuda.flip(cuframe1,1)
        flow_gpu = cuda_farn.calc(cuframe0,cuframe1,None)
    
    frame1 = cv2.cvtColor(frame1,cv2.COLOR_BGR2GRAY)
    frame1 = cv2.flip(frame1,1)
    flow = cv2.calcOpticalFlowFarneback(frame0, frame1, None, .5, 5, of_fb_winsize, 3, 5, 1.1, 0)
    
    return (flow,frame1)


# average angle
def get_angles(flow,frame1,mag,ang):
    # Find the mean vector.
    flow[mag<mag_threshold,0]=np.NaN
    flow[mag<mag_threshold,1]=np.NaN
    X0 = np.nanmean(flow[:,xrange0,0])
    Y0 = np.nanmean(flow[:,xrange0,1])
    if len(xrange1)>0:
        X1 = np.nanmean(flow[:,xrange1,0])
        Y1 = np.nanmean(flow[:,xrange1,1])
    else:
        X1 = 0
        Y1 = 0
        
    if np.isnan(X0):
        X0=0
    if np.isnan(X1):
        X1=0
    if np.isnan(Y0):
        Y0=0
    if np.isnan(Y1):
        Y1=0

    mag0 = np.sqrt(X0**2+Y0**2)
    mag1 = np.sqrt(X1**2+Y1**2)

    ang0 = np.arctan2(Y0,X0)
    if ang0<0:
        ang0 = np.mod(ang0,np.pi)+np.pi
    ang1 = np.arctan2(Y1,X1)
    if ang1<0:
        ang1 = np.mod(ang1,np.pi)+np.pi
    ang2 = np.mod(ang0 - ang1,2*np.pi)

    return (mag0,mag1,ang0,ang1,ang2,X0,Y0,X1,Y1)


def visual_modality_draw_circles(x_targets):
    # Circles moving on the x-axis.
    cv2.circle(vis_mod_bg, (int(vis_mod_center[0]+x_targets[0]),int(vis_mod_center[1]-3)), 30, (250,1,50,1), thickness = 4)
    cv2.circle(vis_mod_bg, (int(vis_mod_center[0]+x_targets[1]),int(vis_mod_center[1]+3)), 20, (200,250,1,1), thickness = 4)
    cv2.imshow('Visual Task', vis_mod_bg)

def visual_modality_draw_lines(theta,radius=100.,total_angle_prop=.7):
    # Lines+tilt.
    # theta[0]=total_angle_prop*theta[0]
    cv2.line(vis_mod_bg, (int(vis_mod_center[0]-radius*np.cos(theta[0])),int(vis_mod_center[1]+radius*np.sin(theta[0]))), (int(vis_mod_center[0]+radius*np.cos(theta[0])),int(vis_mod_center[1]-radius*np.sin(theta[0]))), color=(250,111,150,1), thickness = 7)
    cv2.line(vis_mod_bg, (int(vis_mod_center[0]-radius*np.cos(theta[1])),int(vis_mod_center[1]+radius*np.sin(theta[1]))), (int(vis_mod_center[0]+radius*np.cos(theta[1])),int(vis_mod_center[1]-radius*np.sin(theta[1]))), color=(100,  5,  1,1), thickness = 7)
    cv2.imshow('Visual Task', vis_mod_bg)

def visual_feedback(X0,Y0,X1,Y1,mag):
    # Experiment with the scaling and thresholding to map motion b/w 0 and 255.
    mag[mag<mag_threshold]=0
    mag=mag*10
    
    hsv[...,0] = ang*180/np.pi/2
    hsv[...,2] = mag
    
    # commented this out because you want to be able to tell how much movement is detected, and how fast, from the video.
    #hsv[...,2] = cv2.normalize(mag,None,alpha=0,beta=255,norm_type=cv2.NORM_MINMAX)
    bgr = cv2.cvtColor(hsv,cv2.COLOR_HSV2BGR)
    
    cv2.putText(bgr, datetime.datetime.now().strftime("%A %d %B %Y %I:%M:%S%p"),
                (10, bgr.shape[0] - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.35, (0, 0, 255), 1)
    
    cv2.circle(bgr, center, feedback_circle_r, (25,25,25,1), thickness = 1)
    # Either display individual velocity vectors or the relative phase.
    #if REL_PHASE_FEEDBACK == 1:
    #    cv2.line(bgr, center, (int(center[0] + xrel[0] * feedback_circle_r),int(center[1] + yrel[0] * feedback_circle_r)), (200,200,250,1), thickness = 2)
    #else:
    
    cv2.line(bgr, center, (int(center[0] + (X0) * feedback_circle_r),int(center[1] + (Y0) * feedback_circle_r)), (255,0,255,1), thickness = 2)
    cv2.line(bgr, center, (int(center[0] + (X1) * feedback_circle_r),int(center[1] + (Y1) * feedback_circle_r)), (255,255,0,1), thickness = 2)

    cv2.imshow("Camera", frame1)
    cv2.imshow('Dense optic flow', bgr)

##
if __name__ == '__main__':
    prog_start = time.time()

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv,"hd:e:t:l:p:x:a:i:a:v:s:f:z:",["help","duration=","epsilon=","task=","mov_logging=","mov_plotting=","sound_export=","append_to_file_name=","input=","aud_feedback=","vis_feedback=","synth=","performance_feedback=","visual_modality="])
    except getopt.GetoptError:
        print '\n'
        print 'Some options are not typed correctly. Check with chua_dst_sync_0.97.py -h'
        print '\n'
        sys.exit(2)
    
    # Flag defaults.
    DURATION = 10.0
    MOVT_LOGGING = bool(0)
    MOVT_PLOTTING = bool(1)
    SOUND_EXPORT_FLAG = bool(0)
    AUD_FEEDBACK = bool(1)
    task=1
    APPEND_TO_FILENAME = ''
    mouse_status=True
    wii_status=False
    arduino_status=False
    dt_target=.005 # .005
    dt_ceil = dt_target*1.
    srbar=int(1/dt_target)
    INPUT_DEVICE='mouse'
    RAWAXISLIMIT= 10
    MOUSE_GAIN = 2.
    synth_mode = 1
    FADEIN_R = 1
    FADEOUT_L = 0
    FADEOUT_L_DUR = 10
    EPSILON=.0
    FADEIN_R_TIME = 5
    AMPLIFY_OSC2 = 0
    AMPLIFY_OSC_0 = 20
    AMPLIFY_OSC_0_R = 18
    SCORE_FEEDBACK = bool(0)
    VIS_MODALITY = False
    SOUND_FLAG = True
    cpg_mode='null'
    ForceAdded = .00
    tempo_block_duration = np.Inf
    tempo_block_dur_range = [0,0]
    OMEGA_DRIVER = .5*math.pi
    SCALE_KUR_AMP = .8
    
    """
    ε=1.00 works with Chua1, maybe the task becomes too easy.
    ε was .7 in the first study.
    .5 is OK. .3 is getting hard. .2 = double period.
    .1 Seems impossible to do single harmonic freq.
    """
    

    # Cam parameters
    # Video capture
    CAM_VIS_FEEDBACK=False
    cam_status=False
    GPU_FLAG=False
    feedback_circle_r = 200 # The size of the circle in the center of the screen.
    mag_threshold = .1 # This is important. We try to reduce noise by zero-ing out fluctuations in pixel intensity below a certain threshold.
    FAKE_FRAME_COUNTER = 0

    for opt, arg in opts:
        if opt in ("-h","--help"):
            print'chua_dst_sync_0.95.py --duration=<trial duration> --epsilon=<ε, i.e. .7> --task=<{0:100}> --mov_logging=<1 or 0> --mov_plotting=<1 or 0> --sound_export=<1 or 0> --append_to_file_name=<chars> --input=<mouse (default) or wii or arduino or camera> --aud_feedback=<1 or 0>  --vis_feedback=<1 or 0> --performance_feedback=<1 or 0> --visual_modality=<0 OR 1>'
            print 'where the tasks are:\n'
            print 'Movement sonification only, no stimulus:'
            print '\t  0. _ Sonification only, Stimulus OFF.'
            print 'Periodic generator [Kuramoto]:'
            print '\t 10. _ Periodic non-interactive (sine stimulus ON, movement sonification ON, no interaction).'
            print '\t 11. SCT (like #10, but stimulus fades out, SCT).'
            print '\t 12. _ Periodic non-interactive, different params (like #10, but different amp and duration).'
            print '\t 15. _ Periodic interactive (like #10, but phase-coupled).'
            print '\t 17. Periodic non-interactive, tempo-changing (sine stimulus ON with changing tempo).'
            print 'Non-linear dynamically unstable generator [Chua]:'
            print '\t 20. _ Non-periodic non-interactive (like #25, but with no interaction).'
            print '\t 21. Quasi-periodic non-interactive (driven Chua). [Not finished!]'
            print '\t 25. _ Non-periodic interactive (stimulus, movement sonification, and interaction all ON).'
            print '\t 27. Quasi-periodic non-interactive, tempo-changing (like #17 and #21). [Doesn\'t work yet!]'
            print '\t 29. Non-periodic interactive, w/out mov sonification (like #25).'
            print 'Non-linear dynamically unstable generator [Lorenz]:'
            print '\t 30. _ Non-periodic non-interactive (like #20 but with Lorenz).'
            print 'Dynamically unstable generator [Pole on cart]:'
            print '\t 45. Driven cart-pole system. [Not tested much!]'
            print 'Simulations:'
            print '\t100. Simulations, no wiimote necessary.\n'
            sys.exit()
        if opt in ("-d", "--duration"):
            DURATION = float(arg)
        if opt in ("-e", "--epsilon"):
            EPSILON = float(arg)
        if opt in ("-t", "--task"):
            task = float(arg)
        if opt in ("-l", "--mov_logging"):
            MOVT_LOGGING = bool(float(arg))
        if opt in ("-p", "--mov_plotting"):
            MOVT_PLOTTING = bool(float(arg))
        if opt in ("-x", "--sound_export"):
            SOUND_EXPORT_FLAG = bool(float(arg))
        if opt in ("-a", "--append_to_file_name"):
            APPEND_TO_FILENAME = arg
        if opt in ("-i", "--input"):
            INPUT_DEVICE = arg
        if opt in ("-a", "--aud_feedback"):
            AUD_FEEDBACK = bool(float(arg))
        if opt in ("-v", "--vis_feedback"):
            CAM_VIS_FEEDBACK = bool(float(arg))
        if opt in ("-s", "--synth"):
            synth_mode = float(arg)
        if opt in ("-f", "--performance_feedback"):
            SCORE_FEEDBACK = bool(float(arg))
            if SCORE_FEEDBACK:
                MOVT_LOGGING=bool(1)
        if opt in ("-z", "--visual_modality"):
            VIS_MODALITY = bool(float(arg))
            if VIS_MODALITY:
                SOUND_FLAG = False
                
            
    if synth_mode==0:
        Fs_RATIOS = np.array([[1.],[1.]],dtype=float)
    if synth_mode==1:
        Fs_RATIOS = np.array([[1.,2.,3],[1.,2.,3.]],dtype=float)
    if synth_mode==2:
        Fs_RATIOS = np.array([[1.,1*1.33,3.*1.2],[1,1.875,3.75]],dtype=float)
    if synth_mode==3:
        Fs_RATIOS = np.array([[1,1.78,2.9],[1,1.875,3.75]],dtype=float) # "ooo" and "uuu"
    if synth_mode==4:
        Fs_RATIOS = np.array([[1*1.2,1.9*1.33,3.*1.2],[1,1.875,3.75]],dtype=float) # better?
    #print Fs_RATIOS

    if INPUT_DEVICE=='wii':
        wii_status=True
        mouse_status=False
        CAM_VIS_FEEDBACK = bool(0)

    if INPUT_DEVICE=='arduino':
        wii_status=False
        mouse_status=False
        arduino_status=True
        CAM_VIS_FEEDBACK = bool(0)
        srbar = 500

    if INPUT_DEVICE=='mouse':
        #wii_status=False
        mouse_status=True
        m = PyMouse()
        x_dim, y_dim = m.screen_size()
        CAM_VIS_FEEDBACK = bool(0)

    if INPUT_DEVICE=='camera':
        RAWAXISLIMIT= 2
        cam_status=True
        mouse_status=False
        wii_status=False
        cam_gain = 1
        X_integrated = 0
        if VIS_MODALITY:
            cam_gain = cam_gain*-1


    task_num = 1*task
    # Set task=0 to run simulations.
    if task==100:
        task='SIMULATION'
        wii_status=False
        mouse_status=False
        cpg_mode = 'Chua'
        if 1:
            DRIVER_IN_SIMULATION = 'chua'
            EPSILON2=EPSILON
        else:
            DRIVER_IN_SIMULATION = 'sine'
            EPSILON2=0
        
    # Sonify only. No CPG.
    elif task==0:
        task='SONIFY_ONLY'
        FADEIN_R = 0

    # Sine stimulus. Mov sonification but w/out control. This sets \eps=0.
    elif task==10:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .0
        cpg_mode = 'Kuramoto'
        
    # SCT. Sine stimulus. Mov sonification w/out control.
    elif task==11:
        task='CPG_AND_ACC_AND_SONIFY'
        cpg_mode = 'Kuramoto'
        FADEOUT_L = 1
        FADEIN_R = 0

    # Sine stimulus. Mov sonification but w/out control. This sets \eps=0.
    elif task==12:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .0
        cpg_mode = 'Kuramoto'
        OMEGA_DRIVER = .75*math.pi
        SCALE_KUR_AMP = .6
        
    # Kuramoto stimulus, phase-coupled.
    elif task==15:
        task='CPG_AND_ACC_AND_SONIFY'
        cpg_mode = 'Kuramoto'
        EPSILON = .005
        
    # Sine stimulus with perturbation as changing tempo. 
    # Mov sonification but w/out control. This sets \eps=0.
    elif task==17:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .0
        cpg_mode = 'Kuramoto'
        tempo_block_duration = 1
        tempo_block_dur_range = [1,2]

    elif task==20:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .0
        cpg_mode = 'Chua'
    
    # Simplified Chua stimulus. Mov sonification but w/out control. This sets \eps=0.
    elif task==21:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .0
        cpg_mode = 'ChuaDriven' # use eps=.7 for a double period.
        EPSILON_driven = 1.
        OMEGA_DRIVER = 2.*math.pi

    # That's the interactive, unison task. Stimulus sound, control coupling, and mov sonification.
    elif task==25:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .7
        cpg_mode = 'Chua'

    # Like 4, but with perturbation.
    elif task==27:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .0
        cpg_mode = 'ChuaDriven' # use eps=.7 for a double period.
        EPSILON_driven = 1.
        tempo_block_duration = 5
        tempo_block_dur_range = [5,15]

    # Stimulus sound and coupling to the stimulus, but no mov sonification.
    elif task==29:
        task='CPG_AND_ACC'
        EPSILON = .7
        cpg_mode = 'Chua'
        AMPLIFY_OSC_0_R = 0
    
    elif task==30:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .0
        cpg_mode = 'Lorenz'
        RHO = 28
    
    # That's an interactive task. Stimulus sound, control coupling of an unstable cart system, and mov sonification.
    elif task==45:
        task='CPG_AND_ACC_AND_SONIFY'
        EPSILON = .7
        cpg_mode = 'Cart' # (unstable system), 
        
        
    print '\n'
    print 'TRIAL DURATION', '\t\t', DURATION, 's'
    print 'TASK #', '\t\t\t', int(task_num)
    print 'TASK', '\t\t\t', task
    print 'GENERATOR', '\t\t', cpg_mode
    print 'INPUT DEVICE', '\t\t', INPUT_DEVICE
    print 'ε', '\t\t\t', EPSILON
    print 'MOVT LOGGING', '\t\t', MOVT_LOGGING
    print 'MOVT PLOTTING', '\t\t', MOVT_PLOTTING
    print 'SOUND EXPORT_FLAG', '\t', SOUND_EXPORT_FLAG
    print 'ANY SOUND AT ALL','\t',SOUND_FLAG
    print 'AUDITORY FEEDBACK', '\t', AUD_FEEDBACK
    print 'CAM_FEEDBACK', '\t\t', CAM_VIS_FEEDBACK
    print 'PERFORMANCE FEEDBACK', '\t', SCORE_FEEDBACK
    print 'VISUAL MODALITY', '\t', VIS_MODALITY
    print '\n'


    if cam_status:
        DS_SPEEDUP_0=4
    else:
        DS_SPEEDUP_0=4
    if task==4:
        DS_SPEEDUP_0=4
    if VIS_MODALITY:
        DS_SPEEDUP_0=4
    if cpg_mode=='Lorenz':
        DS_SPEEDUP_0 = 1
        DS_SPEEDUP_0 = DS_SPEEDUP_0*(28/RHO*1)
        
    if cpg_mode=='Kuramoto':
        DS_SPEEDUP_0 = DS_SPEEDUP_0*.8
    DS_SPEEDUP = DS_SPEEDUP_0
    
    #OMEGA_DRIVER = 1*math.pi/5*DS_SPEEDUP
    
    # To set the range of stimulus tempo in task 6.
    omega_range = [OMEGA_DRIVER/2,OMEGA_DRIVER*2.]
    
    if task=='SIMULATION':
        OMEGA_DRIVER = OMEGA_DRIVER/2
        DS_SPEEDUP = 4

    if VIS_MODALITY & (cpg_mode=='Kuramoto'):
        MOUSE_GAIN = -5.
    
    THETA = [.00]*2
    if cpg_mode == 'Lorenz':
        rand = np.random.uniform(2,10,3).astype(float)
    else:
        rand = np.random.uniform(.01,.1,3).astype(float)
    XDST_L = [.00]*3
    XDST_L[0]= float(rand[0])
    XDST_L[1]= float(rand[1])
    XDST_L[2]= float(rand[2])

    rand = np.random.uniform(.01,.1,3)
    XDST_R = [.00]*3
    XDST_R[0]= float(rand[0])
    XDST_R[1]= float(rand[1])
    XDST_R[2]= float(rand[2])

    PI = math.pi
    G = 1.00
    
    theta_smoothed = [0.,0.]
    BARTHETA = [.0]*2

    # Some debugging flags and params
    PRINTING_FLAG        = 0
    PRINT_ACC_FLAG       = 0
    RESCALE_WITH_SIGMOID = 0
    
    # Prepare Wiimote
    SAMPLE_COUNTER = 0
    
    OFFSETS = np.zeros((2,3,2))

    comp = os.uname()
    if comp[1]=='pop-os':
        # The white wii controller!
        OFFSETS[:,:,0] = [[-3.,-5.,-1.],[-107.,99.,100.]]
        # The red wii controller!
        OFFSETS[:,:,1] = [[-24.,-27.,-26.],[-122.,70.,71.]]
        #OFFSETS[:,:,0] = [[-25.,-30.,-25.],[-121.,68.,70.]]
    else:
        # The white wii controller!
        OFFSETS[:,:,1] = [[-3.,-5.,-1.],[-107.,99.,100.]]
        # The red wii controller!
        OFFSETS[:,:,0] = [[-24.,-27.,-26.],[-122.,70.,71.]]
        #OFFSETS[:,:,0] = [[-25.,-30.,-25.],[-121.,68.,70.]]
        
    # The arduino.
    OFFSETS_ARD = np.zeros((2,3))
    OFFSETS_ARD[:,:] = [[336.,339.,333.],[277.,393.,288.]]

    MOVT_BUFFER_SIZE = 12
    X_STATE_BUFFER = [.00] * MOVT_BUFFER_SIZE
    Y_STATE_BUFFER = [.00] * MOVT_BUFFER_SIZE
    Z_STATE_BUFFER = [.00] * MOVT_BUFFER_SIZE
    D_ACC_BUFFER = [.00] * MOVT_BUFFER_SIZE
    x2 = y2 = z2 = 0
    x = y = z = 0
    if 1:
        flip_wii_xaxis_sign = -1
    else:
        flip_wii_xaxis_sign = 1


    EFFECTOR = np.zeros(((srbar * (int(DURATION) + 30)),1),dtype=float)
    DISCRUPDATE = np.zeros(((srbar * (int(DURATION) + 30)),1),dtype=float) 
    FORCEADDED = np.zeros(((srbar * (int(DURATION) + 30)),1),dtype=float)
    THETAVEC = np.zeros(((srbar * (int(DURATION) + 30)),1),dtype=float)
    TIME = np.zeros(((srbar * (int(DURATION) + 30)),1),dtype=float)
    PIXELS = np.zeros(((srbar * (int(DURATION) + 30)),2),dtype=float)
    NOTES = np.zeros(((srbar * (int(DURATION) + 30)),2),dtype=float)
    CPG = np.zeros(((srbar * (int(DURATION) + 30)),3),dtype=float)
    ACC = np.zeros(((srbar * (int(DURATION) + 30)),3,2),dtype=float)
    
    # Prepare sound engine
    if SOUND_FLAG:
        RATE = int(48000/8)
        CHUNK = 8. #8.0 or 20 # ms
        CHUNK = int(CHUNK/1000.0*RATE) # samples
        WAVEDATAL = [0.00] * CHUNK
        WAVEDATAR = [0.00] * CHUNK
        WAVEDATAL2 = [0.00] * CHUNK
        WAVEDATAR2 = [0.00] * CHUNK
        BUFFERS_N = 0
        if SOUND_EXPORT_FLAG:
            sound_frames=[]
    
        # The parameters of the synth oscillator
        BETA  = 10.00**3
        GAMMA = 1
        DT = 1.00/RATE
        XL = np.asarray( np.random.uniform(.01,.1,len(Fs_RATIOS[0])) )
        VL = np.asarray( np.random.uniform(.01,.1,len(Fs_RATIOS[0])) )
        XR = np.asarray( np.random.uniform(.01,.1,len(Fs_RATIOS[0])) )
        VR = np.asarray( np.random.uniform(.01,.1,len(Fs_RATIOS[0])) )
        XL2 = 10.00**(-4)
        VL2 = 10.00**(-4)
        XR2 = 10.00**(-4)
        VR2 = 10.00**(-4)
        FREQ_SQUEEZE = 10
    
        if task=='SONIFY_ONLY':
            AMPLIFY_OSC_L = 0
        else:
            AMPLIFY_OSC_L = 1.*AMPLIFY_OSC_0    
        if AUD_FEEDBACK & ~FADEIN_R:
            AMPLIFY_OSC_R = 1.*AMPLIFY_OSC_0_R
        else:
            AMPLIFY_OSC_R = 0
        
        MIDINOTER = MIDINOTEL = (np.asscalar( np.random.uniform(.01,1,1) ) * 20 + 40)
        NOTES[0,:] = (MIDINOTEL,MIDINOTER)
        FREQL1 = 2.0**((MIDINOTEL - 69) / 12) * 440
        FREQL0 = 1*FREQL1
        FREQL12 = FREQL1/FREQ_SQUEEZE
        FREQL02 = 1*FREQL12    
        FREQR1 = 2.0**((MIDINOTER - 69) / 12) * 440
        FREQR0 = 1*FREQR1
        FREQR12 = FREQR1/FREQ_SQUEEZE
        FREQR02 = 1*FREQR12
        FREQVECL = [FREQL1] * CHUNK
        FREQVECR = [FREQR1] * CHUNK
        FREQVECL2 = [FREQL12] * CHUNK
        FREQVECR2 = [FREQR12] * CHUNK
        
    # Still testing.
    if VIS_MODALITY==True:
        import cv2
        if cpg_mode=='Kuramoto':
            window_size = (800,800)
        else:
            window_size = (1200,400)
        vis_mod_bg = np.zeros((window_size[1],window_size[0],1),dtype='uint8')
        vis_mod_center=(np.int(np.round(window_size[0]/2)),np.int(np.round(window_size[1]/2)))
    

    if cam_status:
        import cv2
        import datetime
        cap = cv2.VideoCapture(0)
        flow,frame0 = first_camera_frame(cap)
        
        # Why download the gpu-ed flow into a var already-processed w/ cpu? Replace w/ if/else?
        if GPU_FLAG:
            flow_gpu.download(flow)
            
        ABS_FRAME_DIFF = np.zeros_like(TIME)
        XY_ABS = np.zeros_like(TIME)
        XY_MAG = np.zeros_like(TIME)
        FAKE_FRAME_COUNTER = 0
        XYmag = np.zeros((len(TIME),3))
        XYang = np.zeros((len(TIME),3))


    if arduino_status:
        import serial
        print 'Important! Using the Arduino for input assumes a specific program has been uploaded.'
        no_serial=True
        if no_serial:
            try:
                ser = serial.Serial('/dev/ttyACM0',baudrate=115200)
                no_serial=False
            except:
                print '\n!!! No Arduino found on /dev/ttyACM0\n'
        if no_serial:
            try:
                ser = serial.Serial('/dev/ttyACM1',baudrate=115200)
                no_serial=False
            except:
                print '\n!!! No Arduino found on /dev/ttyACM0\n'
        ser.close()
        ser.open()
        print 'Print some sample data from the Arduino board:'
        for i in range(0,100):
            data = ser.readline()
        print data

        
    if wii_status:
        import xwiimote
        
        # display a constant
        print "=== " + xwiimote.NAME_CORE + " ==="
    
        # List and remember wiimotes
        try:
            monitor = xwiimote.monitor(True, True)
            # print "Monitor fd", monitor.get_fd(False)
            wiis = list()
            ent = monitor.poll()
            wiis.append(ent)
            while ent is not None:
                print "Found device: " + ent
                ent = monitor.poll()
                wiis.append(ent)
        except SystemError as e:
                print "ooops, cannot create monitor (", e, ")"

        # Continue only if there are wiimotes
        if len(wiis) < 1:
            print "Less than one wiimote to read!"
            exit(0)
        
        # Create two new ifaces
        try:
            devs = list()
            for d in range(0,len(wiis)-1):
                dev = xwiimote.iface(wiis[d])
                devs.append(dev)
        except IOError as e:
            print "ooops,", e
            exit(1)
        
        # Display some information and open the iface
        fds = list()
        for d in range(0,len(devs)):
            f = display_some_info_and_open_iface(devs[d],d+1)
            fds.append(f)
        
        # Read some values
        pwiimote = poll()
        for d in [0,len(devs)-1]:
            pwiimote.register(devs[d].get_fd(), POLLIN)
            
        revt = xwiimote.event()

    # Start the sound engine
    if SOUND_FLAG:
        p = pyaudio.PyAudio()
        devinfo = p.get_device_info_by_index(0)
        NUM_CHANS=devinfo['maxOutputChannels']
        if 0:
            print_sound_device_info(p)
        # print (NUM_CHANS), "channels found."
        if NUM_CHANS<4:
            NUM_CHANS=int(2)
            stream = p.open(format = pyaudio.paFloat32,
                        channels = NUM_CHANS,
                        rate = RATE,
                        output = True,
                        frames_per_buffer=CHUNK,
                        stream_callback = callback2)
        else:
            NUM_CHANS=int(4)
            stream = p.open(format = pyaudio.paFloat32,
                        channels = NUM_CHANS,
                        rate = RATE,
                        output = True,
                        frames_per_buffer=CHUNK,
                        stream_callback = callback4)
        
        stream.start_stream()
    ##    
        
    print '\n'
    print "Start!"
    if SOUND_FLAG:
        if PRINTING_FLAG:
            print "MIDINOTEL:", MIDINOTEL, " @ ", FREQL1, "Hz", "; MIDINOTER:", MIDINOTER, " @ ", FREQR1, "Hz"
    
    start_time      = time.time()
    note_start      = start_time
    last_peak_time  = start_time
    idle_start_time = start_time
    last_tempo_change_time = time.time()
    # Start the trial
    try:
        while True:
            update_states = False
            
            if VIS_MODALITY:
                if SAMPLE_COUNTER > 1:
                    vis_mod_bg = np.multiply(vis_mod_bg,0)
                    if cpg_mode=='Kuramoto':
                        visual_modality_draw_lines(BARTHETA)
                    else:
                        PIXELS[SAMPLE_COUNTER-1,:] = map_x_to_pixel([XDST_L[1]-.5,EFFECTOR[SAMPLE_COUNTER - 1]])
                        visual_modality_draw_circles(PIXELS[SAMPLE_COUNTER-1,:])
                    if cv2.waitKey(1) & 0xFF == ord('q'):
                        break

            if wii_status:
                polls = pwiimote.poll()
                for fd, evt in polls:
                    if fd == fds[0]:
                        devs[0].dispatch(revt)
                        if revt.type == xwiimote.EVENT_ACCEL:
                            x,y,z = revt.get_abs(0)
                            (x,y,z) = rescale_acc((x,y,z))
                            if cpg_mode=='Kuramoto':
                                x = x
                            else:
                                x = flip_wii_xaxis_sign*x
                                # Shift so that the "0" is tilted to the left.
                                x = x+.5
                            update_states = True
                    if len(devs)>1:
                        if fd == fds[1]:
                            devs[1].dispatch(revt)
                            if revt.type == xwiimote.EVENT_ACCEL:
                                x2,y2,z2 = revt.get_abs(0)
                                (x2,y2,z2) = rescale_acc((x2,y2,z2),1)
                                if cpg_mode=='Kuramoto':
                                    x2 = x2
                                else:
                                    x2 = flip_wii_xaxis_sign*x2
                                    # Shift so that the "0" is tilted to the left.
                                    x2 = x2+.5
            
            if arduino_status:
                #if (time.time()-TIME[SAMPLE_COUNTER-1])>dt_target:
                data = ser.readline()
                data=data.split(',')
                if len(data)==4:
                    try:
                        x = int(data[0])
                        y = int(data[1])
                        z = int(data[2])
                        x,y,z = rescale_ard_acc((x,y,z))
                        update_states = True
                    except:
                        print 'Failed to read from Arduino ' + float(np.random.uniform(.01,.1,1))
                        update_states = False

            if mouse_status:
                if (time.time()-TIME[SAMPLE_COUNTER-1])>dt_target:
                    update_states = True
        
                    xy=m.position()
                    if x_dim<2000:
                        xypos=np.array((np.divide(np.ndarray.astype(np.array(xy),'float'),[x_dim,y_dim])-.5)*2)
                    else:
                        xypos=np.array((np.divide(np.ndarray.astype(np.array(xy),'float'),[x_dim/2,y_dim])-.5)*2)
                        
                    # Shift and rescale cursor to be roughly centered at the center of the screen.
                    x=(xypos[0]+.1)*MOUSE_GAIN
                    y=xypos[1]*MOUSE_GAIN
                    if abs(x)>RAWAXISLIMIT:
                        x=x/abs(x)*RAWAXISLIMIT
                    if abs(y)>RAWAXISLIMIT:
                        y=y/abs(y)*RAWAXISLIMIT
            
            
            if task=='SIMULATION':
                if (time.time()-TIME[SAMPLE_COUNTER-1])>dt_target:
                    time.sleep(.0001)
                    update_states = True
                    if SAMPLE_COUNTER > 1:
                        # Scenario 1: Sine wave driver
                        if DRIVER_IN_SIMULATION == 'sine':
                            XDST_R[1] = artificial_SensorGain_sine(OMEGA_DRIVER, TIME[SAMPLE_COUNTER-1] ) #+ np.random.normal(0,.1,1)
                        
                        # Scenario 2: Chua driver
                        if DRIVER_IN_SIMULATION == 'chua':
                            dt = float(time.time()-TIME[SAMPLE_COUNTER-1])
                            
                            EFFECTOR[SAMPLE_COUNTER - 1] = XDST_R[1]
                            XDST_R=Chua1(DS_SPEEDUP*dt,XDST_R,EPSILON2*CPG[SAMPLE_COUNTER-1,1],-1,1,9,5,1) #+ np.random.normal(0,.1,1)
                        x,y,z = XDST_R
            
            
            if cam_status:
                # if the `q` key is pressed, break from the loop
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    break
                
                ret,frame1 = cap.read()
                frame1 = cv2.cvtColor(frame1,cv2.COLOR_BGR2GRAY)
                frame1 = cv2.flip(frame1,1)
                
                if GPU_FLAG:
                    cuframe1 = cv2.cuda_GpuMat(frame1)
            
                # Don't process identical frames, which could happen if the camera is covertly upsampling.
                # Somehow, by coincidence, the fps with flow estimation is just about the 
                # real fps without flow estimation but with skipping identical frames (fake new frames).
                if np.sum(np.abs(frame1-frame0)) == 0:
                    FAKE_FRAME_COUNTER = FAKE_FRAME_COUNTER + 1
                    if np.mod(FAKE_FRAME_COUNTER,100) == 0:
                        print "100 fake frames."
                    continue
                else:
                    update_states = True
            
                
                if GPU_FLAG:
                    cuda_farn.calc(cuframe0,cuframe1,flow_gpu)
                    flow_gpu.download(flow)
                else:
                    flow = cv2.calcOpticalFlowFarneback(frame0, frame1, None, .5, 5, of_fb_winsize, 3, 5, 1.1, 0)
                
                mag, ang = cv2.cartToPolar(flow[...,0], flow[...,1])
                mag0,mag1,ang0,ang1,ang2,x0,y0,x1,y1 = get_angles(flow,frame1,mag,ang)
                
                if CAM_VIS_FEEDBACK:
                    visual_feedback(x0,y0,x1,y1,mag)
                
                XYmag[SAMPLE_COUNTER,0] = mag0
                XYmag[SAMPLE_COUNTER,1] = mag1 
                XYang[SAMPLE_COUNTER,0] = ang0
                XYang[SAMPLE_COUNTER,1] = ang1 
                XYang[SAMPLE_COUNTER,2] = ang2
                
                X_integrated = X_integrated + cam_gain*x0*(time.time() - TIME[SAMPLE_COUNTER - 2])
                if abs(X_integrated)>RAWAXISLIMIT:
                    X_integrated=X_integrated/abs(X_integrated)*RAWAXISLIMIT
                x = X_integrated
                #print x

                frame0 = frame1
                if GPU_FLAG:
                    cuframe0.swap(cuframe1)
                
            #print update_states
            if update_states==True:
                SAMPLE_COUNTER  = SAMPLE_COUNTER + 1
                index_in_buffer = (SAMPLE_COUNTER - 1) % MOVT_BUFFER_SIZE
                ACC[SAMPLE_COUNTER - 1,:,1] = (x2,y2,z2)
                X_STATE_BUFFER[ index_in_buffer ] = ACC[SAMPLE_COUNTER - 1,0,0] = x
                Y_STATE_BUFFER[ index_in_buffer ] = ACC[SAMPLE_COUNTER - 1,1,0] = y
                Z_STATE_BUFFER[ index_in_buffer ] = ACC[SAMPLE_COUNTER - 1,2,0] = z
                TIME[SAMPLE_COUNTER - 1] = time.time() #- start_time

                # Map input devices states to task space variables.
                if task!='SIMULATION':
                    if mouse_status:
                        EFFECTOR[SAMPLE_COUNTER - 1] = X_STATE_BUFFER[ index_in_buffer ]
                        if (cpg_mode=='Kuramoto') & VIS_MODALITY:
                            EFFECTOR[SAMPLE_COUNTER - 1] = convert_to_angle(EFFECTOR - 1.,SAMPLE_COUNTER)
                            #THETA[1] = convert_to_angle(EFFECTOR,SAMPLE_COUNTER)
                        #else:
                            # Even if you force it somehow to work with a mouse, this will never work well with accelerometers.
                            # The noise is nightmare for the phase angle.
                            # Better build a phase-interpolated coupling as in Dotov et al. 2019.
                            #print EFFECTOR[SAMPLE_COUNTER-1], float(convert_to_angle_2(EFFECTOR,SAMPLE_COUNTER))
                            #EFFECTOR[SAMPLE_COUNTER - 1] = convert_to_angle_2(EFFECTOR,SAMPLE_COUNTER)
                    if wii_status:
                        if VIS_MODALITY & (SAMPLE_COUNTER>9):
                            xsmoothed = float(np.mean(ACC[(SAMPLE_COUNTER-10):(SAMPLE_COUNTER-1),0,0]))
                        else:
                            xsmoothed = ACC[SAMPLE_COUNTER-1,0,0]
                        EFFECTOR[SAMPLE_COUNTER - 1] = rescale_x_acc_fun_linear(xsmoothed)
                        if (cpg_mode=='Kuramoto') & VIS_MODALITY:
                            EFFECTOR[SAMPLE_COUNTER - 1] = ((EFFECTOR[SAMPLE_COUNTER - 1]+1)/2*PI)
                        #    EFFECTOR[SAMPLE_COUNTER - 1] = convert_to_angle(EFFECTOR,SAMPLE_COUNTER)
                    if arduino_status:
                        EFFECTOR[SAMPLE_COUNTER - 1] = rescale_x_acc_fun_linear(Y_STATE_BUFFER[index_in_buffer])
                    if cam_status:
                        EFFECTOR[SAMPLE_COUNTER - 1] = X_STATE_BUFFER[index_in_buffer]
                    
                
                if SOUND_FLAG:
                    if task=='SIMULATION':
                        MIDINOTER = map_x_to_note(float(ACC[SAMPLE_COUNTER - 1,1,0]))
                    else:
                        MIDINOTER = map_x_to_note(float(EFFECTOR[SAMPLE_COUNTER - 1]))
                            
                    #print MIDINOTER
                    FREQR1 = midi_key_to_hz(MIDINOTER)
                    FREQR12 = midi_key_to_hz(MIDINOTER)/FREQ_SQUEEZE
                    NOTES  [SAMPLE_COUNTER-1,1] = MIDINOTER

                if task!='SONIFY_ONLY': #task=='CPG_AND_ACC' or task=='CPG_AND_ACC_AND_SONIFY':
                    if SAMPLE_COUNTER > 1:
                        dt=float(time.time()-TIME[SAMPLE_COUNTER-2])
                        dt_ceil = float(min((.05,dt)))
                        if cpg_mode=='Chua':
                            FORCEADDED[SAMPLE_COUNTER - 1]=EPSILON*float(EFFECTOR[SAMPLE_COUNTER-1])
                            XDST_L=Chua1(DS_SPEEDUP*dt_ceil,XDST_L,FORCEADDED[SAMPLE_COUNTER - 1],-1,1,9,5,1)
                                
                        if cpg_mode=='Lorenz':
                            FORCEADDED[SAMPLE_COUNTER - 1]=EPSILON*float(EFFECTOR[SAMPLE_COUNTER-1])
                            XDST_L=Lorenz(DS_SPEEDUP*dt_ceil,XDST_L,FORCEADDED[SAMPLE_COUNTER - 1],RHO)
                            
                        if cpg_mode=='ChuaDriven':
                            FORCEADDED[SAMPLE_COUNTER - 1]=EPSILON_driven*artificial_SensorGain_sine( OMEGA_DRIVER, float(TIME[SAMPLE_COUNTER-1]))
                            XDST_L=Chua1(DS_SPEEDUP*dt_ceil,XDST_L,FORCEADDED[SAMPLE_COUNTER - 1],-1,1,9,5,1)

                        if cpg_mode=='Cart':
                            FORCEADDED[SAMPLE_COUNTER - 1]=float(EFFECTOR[SAMPLE_COUNTER-1])
                            XDST_L[1]=CartPole(dt_ceil,FORCEADDED[SAMPLE_COUNTER - 1],XDST_L[1])
                        
                        if cpg_mode=='Kuramoto':
                            #XDST_L[1]=artificial_SensorGain_sine( OMEGA_DRIVER, TIME[SAMPLE_COUNTER-1] )
                            #FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(np.sin(THETA[1]-THETA[0]))
                            #FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(np.sin(EFFECTOR[SAMPLE_COUNTER-1]-THETA[0]))
                            if VIS_MODALITY:
                                FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(np.cos(THETA[0])-np.cos(BARTHETA[1]))*np.sign(np.sin(THETA[0]))
                                #FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(BARTHETA[1]-BARTHETA[0])*np.sign(EFFECTOR[SAMPLE_COUNTER-1]-EFFECTOR[SAMPLE_COUNTER-2])
                                # This needs either a different task space, not Kuramoto but some kind of spring on a rail, or very smooth sensor.
                                #FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(np.sin(EFFECTOR[SAMPLE_COUNTER-2]-THETA[0]))
                                #FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(np.sin(BARTHETA[1]-THETA[0]-0*PI/2))
                                #FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(BARTHETA[1]-BARTHETA[0])
                                #FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(np.sin(BARTHETA[1]-BARTHETA[0]))
                                #print BARTHETA, FORCEADDED[SAMPLE_COUNTER - 1]
                            else:
                                FORCEADDED[SAMPLE_COUNTER - 1] = EPSILON/2*(MIDINOTER-MIDINOTEL)
                            THETA[0] = (THETA[0] + dt*OMEGA_DRIVER + FORCEADDED[SAMPLE_COUNTER - 1]) % (PI*2)
                            BARTHETA = map_eff_states_to_screen_bars(THETA[0],EFFECTOR[SAMPLE_COUNTER-1])
                            
                            XDST_L[0] = THETA[0]
                            XDST_L[1] = BARTHETA[0]
                            XDST_L[2] = BARTHETA[1]
                            #XDST_L[1] = (np.sin(THETA[0])+1)/2
                            #XDST_L[2] = THETA[1]

                        if (cpg_mode=='Kuramoto') & (tempo_block_duration<np.inf):
                            # Randomly change the stimulus tempo.
                            if 1: # Change at cycle beginning.
                                if SAMPLE_COUNTER>3:
                                    if (np.diff(CPG[(SAMPLE_COUNTER-4):(SAMPLE_COUNTER-2),1])<0) & (np.diff(CPG[(SAMPLE_COUNTER-3):(SAMPLE_COUNTER-1),1])>0):
                                        change_tempo = np.asscalar(np.random.uniform(omega_range[0],omega_range[1],1))
                                        DS_SPEEDUP = DS_SPEEDUP_0*(OMEGA_DRIVER/change_tempo)
                                        OMEGA_DRIVER = change_tempo
                                        #print OMEGA_DRIVER, DS_SPEEDUP
                            else: # Change at given time intervals.
                                if (time.time() - last_tempo_change_time)>tempo_block_duration:
                                    change_tempo = np.asscalar(np.random.uniform(omega_range[0],omega_range[1],1))
                                    DS_SPEEDUP = DS_SPEEDUP_0*(OMEGA_DRIVER/change_tempo)
                                    OMEGA_DRIVER = change_tempo
                                    print OMEGA_DRIVER, DS_SPEEDUP
                                    #print(OMEGA_DRIVER )
                                    last_tempo_change_time = time.time()
                                    # Uniform:
                                    #tempo_block_duration = np.asscalar(np.random.uniform(tempo_block_dur_range[0],tempo_block_dur_range[1],1))
                                    # Add Gaussian to make Brownian
                                    grand = np.random.normal(scale=.2)
                                    tempo_block_duration += grand
                                    if tempo_block_duration < .5:
                                        tempo_block_duration -= grand
                                
                        if cpg_mode=='Lorenz':
                            if np.linalg.norm(XDST_L)>np.Inf:
                                print str(XDST_L)
                                print 'Numerical singularity!'
                                #if SAMPLE_COUNTER>20:
                                #    XDST_L[0] = float(CPG[SAMPLE_COUNTER-20,0])
                                #    XDST_L[1] = float(CPG[SAMPLE_COUNTER-20,1])
                                #    XDST_L[2] = float(CPG[SAMPLE_COUNTER-20,2])
                                #else:
                                XDST_L[0] = float(rand[0])
                                XDST_L[1] = float(rand[1])
                                XDST_L[2] = float(rand[2])
                                
                        if cpg_mode=='Chua':
                            if np.linalg.norm(XDST_L)>2:
                                print str(XDST_L)
                                print 'Numerical singularity!'
                                if SAMPLE_COUNTER>20:
                                    XDST_L[0] = float(CPG[SAMPLE_COUNTER-20,0])
                                    XDST_L[1] = float(CPG[SAMPLE_COUNTER-20,1])
                                    XDST_L[2] = float(CPG[SAMPLE_COUNTER-20,2])
                                else:
                                    XDST_L[0] = float(rand[0])
                                    XDST_L[1] = float(rand[1])
                                    XDST_L[2] = float(rand[2])
                            
                        
                        CPG[SAMPLE_COUNTER-1,:] = XDST_L
                        if SOUND_FLAG:
                            if cpg_mode=='Lorenz':
                                MIDINOTEL = map_x_to_note(XDST_L[2]/40)
                            elif cpg_mode=='Kuramoto':
                                MIDINOTEL = map_x_to_note((np.sin(XDST_L[0])+1)/2*SCALE_KUR_AMP)
                            else:
                                MIDINOTEL = map_x_to_note(XDST_L[1])
                                
                            # Consolidate into the vars that control the sound, and log
                            #FORCEADDED[SAMPLE_COUNTER - 1] = ForceAdded
                            FREQL1 = midi_key_to_hz(MIDINOTEL)
                            FREQL12 = FREQL1/FREQ_SQUEEZE
                            NOTES[SAMPLE_COUNTER-1,0] = MIDINOTEL

                if PRINT_ACC_FLAG:
                    print "ACC", "%8.2f" % x, "%8.2f" % y, "%8.2f" % z, "%8.2f" % D_ACC_BUFFER[ index_in_buffer ], "%8.2f" % CPG[ index_in_buffer,0], "%8.2f" % CPG[ index_in_buffer,1], "%8.2f" % CPG[index_in_buffer,2]
            
            
            # Fade-in R
            if SOUND_FLAG:
                if FADEIN_R & AUD_FEEDBACK & ((time.time() - start_time) > FADEIN_R_TIME):
                    AMPLIFY_OSC_R = AMPLIFY_OSC_R + .005*(AMPLIFY_OSC_0_R-AMPLIFY_OSC_R)
                    if AMPLIFY_OSC_R == AMPLIFY_OSC_0_R:
                        FADEIN_R = 0
            
            if SOUND_FLAG:
                if FADEOUT_L & AUD_FEEDBACK & ((time.time() - start_time) > FADEOUT_L_DUR):
                    AMPLIFY_OSC_L = AMPLIFY_OSC_L + .005*(0-AMPLIFY_OSC_L)
                    if AMPLIFY_OSC_L == 0:
                        FADEOUT_L = 0

            if (time.time() - start_time) > DURATION:
                break
            
    except KeyboardInterrupt:
        print "exiting..."

    # Start cleaning up.
    if wii_status:
        for d in range(0,len(devs)):
            pwiimote.unregister(devs[d].get_fd())
            devs[d].set_led(1, True)
            devs[d].set_led(2, False)
            devs[d].set_led(3, False)
            devs[d].set_led(4, (dev.get_battery() < 50))
    
    # Stop the sound stream
    if SOUND_FLAG:
        stream.stop_stream()
        stream.close()
        
    # cleanup the camera and close any open windows
    if cam_status:
        cv2.destroyAllWindows()
        cap.release()
        time.sleep(.1)

    if VIS_MODALITY:
        cv2.destroyAllWindows()
        time.sleep(.1)

    if SOUND_FLAG:
        p.terminate()
    
    trial_duration = time.time() - start_time
    FPS = (SAMPLE_COUNTER / trial_duration)
    print "... Done!"
    print '\n'
    print "Trial + prep time ", '\t', "%.3f" % (time.time() - prog_start), "s"
    print "Trial duration ", '\t', "%.3f" % trial_duration, "s"
    print "Samples collected", '\t', SAMPLE_COUNTER
    print "<SR>", '\t\t\t', "%.2f" % FPS, "Hz"
    
    TIME = TIME[1:SAMPLE_COUNTER-1]
    TIME = TIME-min(TIME)
    ACC = ACC[1:SAMPLE_COUNTER-1,:,:]
    EFFECTOR = EFFECTOR[1:SAMPLE_COUNTER-1]
    CPG = CPG[1:SAMPLE_COUNTER-1,:]
    NOTES = NOTES[1:SAMPLE_COUNTER-1,:]
    PIXELS= PIXELS[1:SAMPLE_COUNTER-1,:]
    DISCRUPDATE = DISCRUPDATE[1:SAMPLE_COUNTER-1]
    FORCEADDED = FORCEADDED[1:SAMPLE_COUNTER-1]
    
    # Export the movement
    if MOVT_LOGGING:
        EPSILON = int(EPSILON*1e2)
        log_file_name = 'trial_log-' + time.strftime("%y%m%d-%H%M%S") + '_task' + "%03.0f" % task_num + '_aud' + "%1d" % SOUND_FLAG + '_vis' + "%1d" % VIS_MODALITY + '_eps' + "%03d" % EPSILON + APPEND_TO_FILENAME
        f = open(log_file_name,'w')
        f.write("%10s," % 'Time'+"%10s," % 'Acc1X'+"%10s," % 'Acc1Y'+"%10s," % 'Acc1Z' + "%10s," % 'AbsAcc'+"%8s," % 'EFFECTOR'+"%10s," % 'CPG'+"%10s," % 'CPGY'+"%10s," % 'CPGZ'+"%10s," % 'Force'+"%8s," % 'NoteL'+"%8s," % 'NoteR'+"%10s," % 'Acc2X'+"%10s," % 'Acc2Y' + "%10s," % 'Acc2Z' + "%10s," % 'XPixStim' + "%10s," % 'XPixPart' + "%12s\n" % 'DiscrUpdate')
        for n in range(0, len(ACC), 1):
            f.write("%10.3f," % TIME[n]+"%10.4f," % ACC[n,0,0]+"%10.4f," % ACC[n,1,0] + "%10.4f," % ACC[n,2,0] + "%10.4f," % 0 +"%8.3f," % EFFECTOR[n]+"%10.3f," % CPG[n,0]+"%10.3f," % CPG[n,1]+"%10.3f," % CPG[n,2]+"%10.3f," % FORCEADDED[n]+"%8.2f," % NOTES[n,0]+"%8.2f," % NOTES[n,1] +"%10.4f," % ACC[n,0,1] + "%10.4f," % ACC[n,1,1] + "%10.4f," % ACC[n,2,1] + "%10.4f," % PIXELS[n,0] + "%10.4f," % PIXELS[n,1] + "%12.2f\n" % DISCRUPDATE[n])
        f.close()
    
    if SCORE_FEEDBACK:
        if MOVT_LOGGING:
            if VIS_MODALITY:
                compute_and_return_performance_feedback(log_file_name,'XPixStim','XPixPart')
            else:
                compute_and_return_performance_feedback(log_file_name,'NoteL','NoteR')
        else:
            print 'You need to be saving (--mov_logging=1) the trials to give performance feedback.'
    
    # Plot
    if MOVT_PLOTTING:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10,10))

        if SOUND_FLAG:
            plt.plot(TIME,NOTES[:,0],'-',label='Note Stim')
            plt.plot(TIME,NOTES[:,1],'-',label='Note Participant')
            plt.ylabel('MIDI Note')
            plt.xlabel('Time, s')
            plt.legend(loc='upper right', shadow=False, fontsize='small')    
            plt.show()
            
        if VIS_MODALITY:
            plt.plot(TIME,PIXELS[:,0],'-',label='X Stim')
            plt.plot(TIME,PIXELS[:,1],'-',label='X Participant')
            plt.ylabel('X')
            plt.xlabel('Time, s')
            plt.legend(loc='upper right', shadow=False, fontsize='small')    
            plt.show()
            
        if wii_status:
            if len(devs)>1:
                plt.plot(TIME,ACC[:,0,0],'-',label='Wii-1, X')
                plt.plot(TIME,ACC[:,0,1],'-',label='Wii-2, X')
                plt.xlabel('Time, s')
                plt.ylabel('MIDI Note')
                plt.legend(loc='upper right', shadow=False, fontsize='small')    
                plt.show()

        # Visualize the raw data of the stimulus and the EFFECTORs.
        if True:
            plt.plot(TIME,ACC[:,0,0],'-',label='Participant X')
            plt.plot(TIME,ACC[:,1,0],'-',label='Participant Y')
            plt.plot(TIME,ACC[:,2,0],'-',label='Participant Z')
            plt.plot(TIME,EFFECTOR,'-',label='Motor Output in task space (inclination, x-dim, etc.)')
            plt.plot(TIME,np.multiply(FORCEADDED,10),'-',label='The coupling function from participant to stimulus')
            plt.plot(TIME,CPG[:,0],'-',label='Stim X')
            plt.plot(TIME,CPG[:,1],'-',label='Stim Y')
            plt.plot(TIME,CPG[:,2],'-',label='Stim Z')
            plt.xlabel('Time, s')
            plt.ylabel('X')
            plt.legend(loc='upper right', shadow=False, fontsize='small')    
            plt.show()            
    
        if False:
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure(figsize=(10,10))
            ax = fig.gca(projection='3d')
            ax.plot(CPG[300:,0],CPG[300:,1],CPG[300:,2],'-k')
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            ax.set_title("Driven Chua Attractor")
            plt.show()

    # Export the sound
    if SOUND_EXPORT_FLAG:
        wavefile=wave.open('sound_output' + time.strftime("%m-%d-%y-%H-%M-%S") + '.wav','w')
        wavefile.setframerate(RATE)
        wavefile.setsampwidth(2)
        wavefile.setnchannels(2)
        wavefile.setnframes(0)
        wavefile.writeframesraw(np.asarray(np.uint16(np.multiply(sound_frames,255**2 ))).tostring())
        wavefile.close()